#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <unistd.h>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_filemanager = new FileManager;
    m_centerline = new Centerline;
    m_colon = new Colon;
    m_rendermanager = new RenderManager;

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->AddRenderer(m_rendermanager->GetRender());

    m_lightdialog.setscene(m_rendermanager);
    m_lightdialog.setwindow(this);
}

MainWindow::~MainWindow()
{
    delete ui;
    delete m_filemanager;
    delete m_centerline;
    delete m_colon;
    delete m_rendermanager;
}
void MainWindow::addlight()
{
    m_rendermanager->GetRender()->RemoveAllLights();
    m_rendermanager->addlight();
}

// load colon surface
void MainWindow::on_actionNew_file_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"), "",
                tr("3Dmodels (*.stl *.vtp *.ply *.obj *.off)"));
    if(filePath.isEmpty()) return;

    m_filemanager->LoadNewFile(filePath);
    m_colon->Object::SetInput(m_filemanager->getfile());
    m_rendermanager->GetRender()->SetBackground(0.5,0.5, 1);
    addlight();
    //m_filemanager->SaveFile(m_colon->GetOutput(), "TexturedColon.vtp");
    m_rendermanager->renderModel(m_colon->GetActor());

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

// centerline and colon deformation are done in this function
void MainWindow::on_actionLoad_Centerline_triggered()
{
    //m_centerline->ConnectTwoContoursTest(m_rendermanager, m_filemanager);

    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"),"",
                tr("Centerline File (*.vtp)"));
    if(filePath.isEmpty()) return;
    m_filemanager->LoadNewFile(filePath);
    m_centerline->Object::SetInput(m_filemanager->getfile());

    // Uniform Sampling
    m_centerline->UniformSample(1000);

    // Gaussian Smoothing
    //m_centerline->SmoothCenterline(3);
    //m_filemanager->SaveFile(m_centerline->GetOutput(), "SmoothedCenterline.vtp");

    // Centerline-Driven Colon Deformation
    /*
    vtkSmartPointer<vtkPolyData> newColonPoly = vtkSmartPointer<vtkPolyData>::New();
    newColonPoly = m_centerline->EliminateTorsion(m_rendermanager, m_colon->GetOutput(), m_filemanager);
    m_filemanager->SaveFile(m_centerline->GetOutput(), "ModifiedCenterline.vtp");
    */
    //m_colon->SetPoint(newColonPoly->GetPoints());


    m_rendermanager->renderModel(m_centerline->GetActor());

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();

}

// Give two points on colon surface, find the geodesic path between them.
// This is used to find the slitting line on the surface
vtkSmartPointer<vtkPolyData> MainWindow::GeodesicPath(vtkSmartPointer<vtkPolyData> points)
{
    vtkSmartPointer<vtkPolyData> path = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> path_temp = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();

    dijkstra->SetInputData(m_colon->GetOutput());

    vtkSmartPointer<vtkPolyData> colon_temp = vtkSmartPointer<vtkPolyData>::New();
    colon_temp = m_colon->GetOutput();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(m_colon->GetOutput());
    pointLocator->BuildLocator();
    for(vtkIdType i = 0; i<points->GetNumberOfPoints()-1; i++)
    {
        vtkIdType StartVertex, EndVertex;
        double closestPoint1[3], closestPoint2[3];
        double point1[3];
        points->GetPoint(i, point1);
        double point2[3];
        points->GetPoint(i+1, point2);

        StartVertex = pointLocator->FindClosestPoint(point1);
        EndVertex = pointLocator->FindClosestPoint(point2);
        colon_temp->GetPoint(StartVertex, closestPoint1);
        colon_temp->GetPoint(EndVertex, closestPoint2);

        //std::cout<<point1[0]<<" "<<point1[1]<<" "<<point1[2]<<std::endl;
        //std::cout<<point2[0]<<" "<<point2[1]<<" "<<point2[2]<<std::endl;
        //std::cout<<closestPoint1[0]<<" "<<closestPoint1[1]<<" "<<closestPoint1[2]<<std::endl;
        //std::cout<<closestPoint2[0]<<" "<<closestPoint2[1]<<" "<<closestPoint2[2]<<std::endl;
        //std::cout<<"from  "<<StartVertex<<" to  "<<EndVertex<<std::endl;

        dijkstra->SetStartVertex(StartVertex);
        dijkstra->SetEndVertex(EndVertex);
        dijkstra->Update();
        path_temp->ShallowCopy(dijkstra->GetOutput());

        appendFilter->RemoveAllInputs();
        appendFilter->AddInputData(path);
        appendFilter->AddInputData(path_temp);
        appendFilter->Update();

        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();

        path->DeepCopy(cleanFilter->GetOutput());
    }
    return path;
}

// Upsample the path by factor 10. Currently unused.
vtkSmartPointer<vtkPolyData> MainWindow::Upsampling(vtkSmartPointer<vtkPolyData> path)
{
    vtkSmartPointer<vtkCardinalSpline> spline = vtkSmartPointer<vtkCardinalSpline>::New();
    vtkSmartPointer<vtkSplineFilter> splineFilter = vtkSmartPointer<vtkSplineFilter>::New();
    splineFilter->SetInputData(path);
    splineFilter->SetNumberOfSubdivisions(path->GetNumberOfPoints() * 10);
    splineFilter->SetSpline(spline);
    splineFilter->Update();
    return splineFilter->GetOutput();
}

void MainWindow::on_action_Computing_triggered()
{ 
    vtkSmartPointer<vtkPoints> points_1 = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> points_c = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> points_2 = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> polypoints_1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> polypoints_c = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> polypoints_2 = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    double closestPoint[3], oppositeClosestPoint[3], currentCenterPoint[3], oppositeFarPoint[3];
    double closestDist;
    vtkIdType cellID;
    int subID;

    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter_1 = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter_2 = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter_C = vtkSmartPointer<vtkVertexGlyphFilter>::New();

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetInputData(m_colon->GetOutput());

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    double normalDirection[3];
    double *p1, *p2, *intersectionPoint;
    p1 = new double [3];
    p2 = new double [3];
    intersectionPoint = new double [3];
    double intersectionT;

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(cutter->GetOutputPort());
    connectivityFilter->SetExtractionModeToClosestPointRegion();

    vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cutterMapper->SetInputConnection(connectivityFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
    planeActor->GetProperty()->SetColor(2, 0.5, 0);
    planeActor->GetProperty()->SetLineWidth(4);
    planeActor->SetMapper(cutterMapper);

    m_rendermanager->renderModel(planeActor);

    QVTKWidget* widget = this->findChild<QVTKWidget*>("qvtk");
    for(vtkIdType i=0; i<m_centerline->GetNumberOfPoints(); i += 5)
    {
        connectivityFilter->SetClosestPoint(m_centerline->GetOutput()->GetPoint(i));
        plane = m_centerline->GetVerticalPlane(i, normalDirection);
        cutter->SetCutFunction(plane);
        cutter->Update();
        connectivityFilter->Update();
        cutterMapper->Update();

        cellLocator->SetDataSet(connectivityFilter->GetOutput());
        cellLocator->BuildLocator();

        m_centerline->GetOutput()->GetPoint(i, currentCenterPoint);
        if(i == 0)
        {
            cellLocator->FindClosestPoint(currentCenterPoint,closestPoint,cellID,subID,closestDist);
            points_1->InsertNextPoint(closestPoint);
            p1[0] = closestPoint[0];
            p1[1] = closestPoint[1];
            p1[2] = closestPoint[2];
            vtkMath::Add(p1, normalDirection, p2);     
        }
        else
        {
            plane->IntersectWithLine(p1, p2, intersectionT, intersectionPoint);
            cellLocator->FindClosestPoint(intersectionPoint, closestPoint, cellID, subID, closestDist);
            points_1->InsertNextPoint(closestPoint);
            p1[0] = closestPoint[0];
            p1[1] = closestPoint[1];
            p1[2] = closestPoint[2];
            vtkMath::Add(p1, normalDirection, p2);         
        }

        oppositeFarPoint[0] = closestPoint[0] + 10 * (currentCenterPoint[0] - closestPoint[0]);
        oppositeFarPoint[1] = closestPoint[1] + 10 * (currentCenterPoint[1] - closestPoint[1]);
        oppositeFarPoint[2] = closestPoint[2] + 10 * (currentCenterPoint[2] - closestPoint[2]);

        cellLocator->FindClosestPoint(oppositeFarPoint, oppositeClosestPoint, cellID, subID, closestDist);
        points_2->InsertNextPoint(oppositeClosestPoint);

        points_c->InsertNextPoint(m_centerline->GetOutput()->GetPoint(i));
        widget->GetRenderWindow()->Render();
    }
    polypoints_1->SetPoints(points_1);
    vertexFilter_1->SetInputData(polypoints_1);
    vertexFilter_1->Update();
    vtkSmartPointer<vtkPolyData> PolyPoints_1 = vtkSmartPointer<vtkPolyData>::New();
    PolyPoints_1 = vertexFilter_1->GetOutput();

    polypoints_2->SetPoints(points_2);
    vertexFilter_2->SetInputData(polypoints_2);
    vertexFilter_2->Update();
    vtkSmartPointer<vtkPolyData> PolyPoints_2 = vtkSmartPointer<vtkPolyData>::New();
    PolyPoints_2 = vertexFilter_2->GetOutput();

    polypoints_c->SetPoints(points_c);
    vertexFilter_C->SetInputData(polypoints_c);
    vertexFilter_C->Update();
    vtkSmartPointer<vtkPolyData> PolyPoints_C = vtkSmartPointer<vtkPolyData>::New();
    PolyPoints_C = vertexFilter_C->GetOutput();


    vtkSmartPointer<vtkPolyDataMapper> pointsMapper_1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> pointsActor_1 = vtkSmartPointer<vtkActor>::New();
    pointsActor_1->SetMapper(pointsMapper_1);
    pointsActor_1->GetProperty()->SetPointSize(4);
    pointsActor_1->GetProperty()->SetOpacity(1);
    pointsActor_1->GetProperty()->SetColor(0,0,1);
    pointsMapper_1->SetInputConnection(vertexFilter_1->GetOutputPort());
    pointsMapper_1->Update();
    m_rendermanager->renderModel(pointsActor_1);
    widget->GetRenderWindow()->Render();

    vtkSmartPointer<vtkPolyDataMapper> pointsMapper_2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> pointsActor_2 = vtkSmartPointer<vtkActor>::New();
    pointsActor_2->SetMapper(pointsMapper_2);
    pointsActor_2->GetProperty()->SetPointSize(4);
    pointsActor_2->GetProperty()->SetOpacity(1);
    pointsActor_2->GetProperty()->SetColor(0,1,0);
    pointsMapper_2->SetInputConnection(vertexFilter_2->GetOutputPort());
    pointsMapper_2->Update();
    m_rendermanager->renderModel(pointsActor_2);
    widget->GetRenderWindow()->Render();

    vtkSmartPointer<vtkPolyData> path_1 = vtkSmartPointer<vtkPolyData>::New();
    path_1 = GeodesicPath(PolyPoints_1);
    vtkSmartPointer<vtkPolyDataMapper> pathMapper_1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    pathMapper_1->SetInputData(path_1);
    vtkSmartPointer<vtkActor> pathActor_1 = vtkSmartPointer<vtkActor>::New();
    pathActor_1->SetMapper(pathMapper_1);
    pathActor_1->GetProperty()->SetColor(1,0,0);
    pathActor_1->GetProperty()->SetLineWidth(5);
    std::cout<<"path_1: "<<path_1->GetNumberOfPoints()<<endl;
    m_rendermanager->renderModel(pathActor_1);
    widget->GetRenderWindow()->Render();

    vtkSmartPointer<vtkPolyData> path_2 = vtkSmartPointer<vtkPolyData>::New();
    path_2 = GeodesicPath(PolyPoints_2);
    vtkSmartPointer<vtkPolyDataMapper> pathMapper_2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    pathMapper_2->SetInputData(path_2);
    vtkSmartPointer<vtkActor> pathActor_2 = vtkSmartPointer<vtkActor>::New();
    pathActor_2->SetMapper(pathMapper_2);
    pathActor_2->GetProperty()->SetColor(0,0,1);
    pathActor_2->GetProperty()->SetLineWidth(5);
    std::cout<<"path_2: "<<path_2->GetNumberOfPoints()<<endl;
    m_rendermanager->renderModel(pathActor_2);
    widget->GetRenderWindow()->Render();

    vtkSmartPointer<vtkPolyData> colon_temp = vtkSmartPointer<vtkPolyData>::New();
    colon_temp = m_colon->GetOutput();

    // upsample the path_1 and path_2
    /*
    path_1->DeepCopy(Upsampling(path_1));
    std::cout<<"path_1: "<<path_1->GetNumberOfPoints()<<endl;
    path_2->DeepCopy(Upsampling(path_2));
    std::cout<<"path_2: "<<path_2->GetNumberOfPoints()<<endl;
    */

    //select the points along the path out, then inverse, get the mesh with out the path area
    /*
    vtkSmartPointer<vtkIdTypeArray> idsPath = vtkSmartPointer<vtkIdTypeArray>::New();

    vtkSmartPointer<vtkPointLocator> locatorPath = vtkSmartPointer<vtkPointLocator>::New();
    locatorPath->SetDataSet(m_colon->GetOutput());
    locatorPath->BuildLocator();


    double firstTangent[3], lastTangent[3], tempdir1[3], tempdir2[3], firstCPoint[3], lastCPoint[3];
    double distanceThreshold = 2500;
    m_centerline->GetOutput()->GetPoint(50, firstCPoint);
    m_centerline->GetOutput()->GetPoint(m_centerline->GetNumberOfPoints()-50, lastCPoint);
    m_centerline->GetVerticalPlane(50, firstTangent);
    m_centerline->GetVerticalPlane(m_centerline->GetNumberOfPoints()-50, lastTangent);

    for(vtkIdType i=0; i<m_colon->GetOutput()->GetNumberOfPoints(); i++)
    {
        double p[3];
        m_colon->GetOutput()->GetPoint(i, p);
        // remove the head and tail of colon
        vtkMath::Subtract(p, firstCPoint, tempdir1);
        if(vtkMath::Dot(tempdir1, firstTangent) < 0 && vtkMath::Distance2BetweenPoints(p, firstCPoint) < distanceThreshold * 2)
            idsPath->InsertNextValue(i);
        vtkMath::Subtract(p, lastCPoint, tempdir2);
        if(vtkMath::Dot(tempdir2, lastTangent) > 0 && vtkMath::Distance2BetweenPoints(p, lastCPoint) < distanceThreshold * 4)
            idsPath->InsertNextValue(i);
    }

    for(vtkIdType i = 0; i< path_1->GetNumberOfPoints(); i++)
    {
        double tp[3];
        vtkIdType tid;
        path_1->GetPoint(i, tp);
        tid = locatorPath->FindClosestPoint(tp);
        idsPath->InsertNextValue(tid);
    }
    for(vtkIdType i = 0; i< path_2->GetNumberOfPoints(); i++)
    {
        double tp[3];
        vtkIdType tid;
        path_2->GetPoint(i, tp);
        tid = locatorPath->FindClosestPoint(tp);
        idsPath->InsertNextValue(tid);
    }

    vtkSmartPointer<vtkSelectionNode> selectPathNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectPathNode->SetFieldType(vtkSelectionNode::POINT);
    selectPathNode->SetContentType(vtkSelectionNode::INDICES);
    selectPathNode->SetSelectionList(idsPath);
    selectPathNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);
    selectPathNode->GetProperties()->Set(vtkSelectionNode::INVERSE(), 1);

    vtkSmartPointer<vtkSelection> selectPath = vtkSmartPointer<vtkSelection>::New();
    selectPath->AddNode(selectPathNode);

    vtkSmartPointer<vtkExtractSelection> extractselectPath = vtkSmartPointer<vtkExtractSelection>::New();
    extractselectPath->SetInputData(0, m_colon->GetOutput());
    extractselectPath->SetInputData(1, selectPath);
    extractselectPath->Update();

    vtkSmartPointer<vtkUnstructuredGrid> selectedSurface = vtkSmartPointer<vtkUnstructuredGrid>::New();
    selectedSurface->ShallowCopy(extractselectPath->GetOutput());

    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(selectedSurface);
    geometryFilter->Update();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilterForHalfColon = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connectivityFilterForHalfColon->SetInputConnection(geometryFilter->GetOutputPort());
    connectivityFilterForHalfColon->SetExtractionModeToLargestRegion();
    connectivityFilterForHalfColon->Update();

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputConnection(connectivityFilterForHalfColon->GetOutputPort());

    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    m_showselectedwindow.show();
    m_showselectedwindow.RenderSelected(selectedActor);
    */

    // selection of half the colon the very first version
    /*
    vtkSmartPointer<vtkCellLocator> cellLocator2 = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator2->SetDataSet(PolyPoints_C);
    cellLocator2->BuildLocator();
    vtkSmartPointer<vtkPolyData> colon_temp = vtkSmartPointer<vtkPolyData>::New();
    colon_temp = m_colon->GetOutput();


    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    double closestPC[3], p[3], pnext[3], plast[3];
    double direction[3], direction_last[3], direction_next[3], directionV[3], directionN[3], directionP[3];
    for(vtkIdType i=0; i<colon_temp->GetNumberOfPoints(); i++)
    {
        colon_temp->GetPoint(i, p);
        cellLocator2->FindClosestPoint(p, closestPC, cellID, subID, closestDist);
        double validationp[3];
        PolyPoints_C->GetPoint(cellID, validationp);
        //std::cout<<i<<"->"<<cellID << ": " << closestPC[0]<<" "<<closestPC[1]<<" "<<closestPC[2]<<" "<< validationp[0]<<" "<<validationp[1]<<" "<<validationp[2]<<" ";

        if(i == 0)
        {
            PolyPoints_C->GetPoint(cellID + 1 , pnext);
            vtkMath::Subtract(pnext, p, direction);
            vtkMath::Normalize(direction);
        }
        else if(i == colon_temp->GetNumberOfPoints()-1)
        {
            PolyPoints_C->GetPoint(i-1 , plast);
            vtkMath::Subtract(p, plast, direction);
            vtkMath::Normalize(direction);
        }
        else
        {
            PolyPoints_C->GetPoint(cellID + 1, pnext);
            PolyPoints_C->GetPoint(cellID - 1, plast);
            vtkMath::Subtract(pnext, p, direction_next);
            vtkMath::Normalize(direction_next);
            vtkMath::Subtract(p, plast, direction_last);
            vtkMath::Normalize(direction_last);
            vtkMath::Add(direction_next, direction_last, direction);
            vtkMath::Normalize(direction);
        }
        vtkMath::Subtract(PolyPoints_1->GetPoint(cellID), closestPC, directionV);
        vtkMath::Normalize(directionV);
        vtkMath::Cross(direction, directionV, directionN);

        vtkMath::Subtract(p, closestPC, directionP);
        if(vtkMath::Dot(directionP, directionN) >=0)
        {
            //std::cout<<"yes"<<std::endl;
            ids->InsertNextValue(i);
        }
        //else
            //std::cout<<"no"<<std::endl;
    }

    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);
    selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, m_colon->GetOutput());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(selected);
    geometryFilter->Update();

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(geometryFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(50);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputConnection(smoothFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    m_showselectedwindow.show();
    m_showselectedwindow.RenderSelected(selectedActor);
    */


    // selection of half of the colon
    vtkSmartPointer<vtkPointLocator> pointLocator_C = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator_C->SetDataSet(m_centerline->GetOutput());
    pointLocator_C->BuildLocator();
    vtkSmartPointer<vtkPointLocator> pointLocator_1 = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator_1->SetDataSet(path_1);
    pointLocator_1->BuildLocator();
    vtkSmartPointer<vtkPointLocator> pointLocator_2 = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator_2->SetDataSet(path_2);
    pointLocator_2->BuildLocator();

    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);
    double closestPC[3], closestP1[3], closestP2[3], p[3], pnext[3], plast[3];
    double direction[3], direction_last[3], direction_next[3], directionV[3], directionN[3], directionP[3];
    vtkIdType Id1, Id2, IdC;


    double firstTangent[3], lastTangent[3], tempdir1[3], tempdir2[3], firstCPoint[3], lastCPoint[3];
    double distanceThreshold = 2500;
    m_centerline->GetOutput()->GetPoint(0, firstCPoint);
    m_centerline->GetOutput()->GetPoint(m_centerline->GetNumberOfPoints()-1, lastCPoint);
    m_centerline->GetVerticalPlane(0, firstTangent);
    m_centerline->GetVerticalPlane(m_centerline->GetNumberOfPoints()-1, lastTangent);

    for(vtkIdType i=0; i<colon_temp->GetNumberOfPoints(); i++)
    {
        colon_temp->GetPoint(i, p);

        // remove the head and tail of colon
        vtkMath::Subtract(p, firstCPoint, tempdir1);
        if(vtkMath::Dot(tempdir1, firstTangent) < 0 && vtkMath::Distance2BetweenPoints(p, firstCPoint) < distanceThreshold)
            continue;
        vtkMath::Subtract(p, lastCPoint, tempdir2);
        if(vtkMath::Dot(tempdir2, lastTangent) > 0 && vtkMath::Distance2BetweenPoints(p, lastCPoint) < distanceThreshold)
            continue;

        IdC = pointLocator_C->FindClosestPoint(p);
        m_centerline->GetOutput()->GetPoint(IdC, closestPC);
        Id1 = pointLocator_1->FindClosestPoint(p);
        path_1->GetPoint(Id1, closestP1);
        Id2 = pointLocator_2->FindClosestPoint(p);
        path_2->GetPoint(Id2, closestP2);

        if(i == 0)
        {
            m_centerline->GetOutput()->GetPoint(IdC + 1 , pnext);
            vtkMath::Subtract(pnext, p, direction);
            vtkMath::Normalize(direction);
        }
        else if(i == colon_temp->GetNumberOfPoints()-1)
        {
            m_centerline->GetOutput()->GetPoint(IdC - 1 , plast);
            vtkMath::Subtract(p, plast, direction);
            vtkMath::Normalize(direction);
        }
        else
        {
            m_centerline->GetOutput()->GetPoint(IdC + 1, pnext);
            m_centerline->GetOutput()->GetPoint(IdC - 1, plast);
            vtkMath::Subtract(pnext, p, direction_next);
            vtkMath::Normalize(direction_next);
            vtkMath::Subtract(p, plast, direction_last);
            vtkMath::Normalize(direction_last);
            vtkMath::Add(direction_next, direction_last, direction);
            vtkMath::Normalize(direction);
        }
        vtkMath::Subtract(closestP1, closestP2, directionV);
        vtkMath::Normalize(directionV);
        vtkMath::Cross(direction, directionV, directionN);

        vtkMath::Subtract(p, closestPC, directionP);
        if(vtkMath::Dot(directionP, directionN) >=0)
            ids->InsertNextValue(i);
    }

    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);
    selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, m_colon->GetOutput());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(selected);
    geometryFilter->Update();

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(geometryFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(50);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputConnection(smoothFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    m_showselectedwindow.show();
    m_showselectedwindow.RenderSelected(selectedActor);


    delete p1;
    delete p2;
    delete intersectionPoint;
    std::cout<<"computing finished"<<std::endl;
}

void MainWindow::on_actionLighting_triggered()
{
    m_lightdialog.exec();
}

void MainWindow::on_action_Deform_Colon_triggered()
{
    /*
    vtkSmartPointer<vtkPolyData> newcenterline = vtkSmartPointer<vtkPolyData>::New();
    newcenterline = m_centerline->EliminateTorsion(m_rendermanager, m_colon->GetOutput(), m_filemanager);


    vtkSmartPointer<vtkVertexGlyphFilter> newVertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    newVertexFilter->SetInputData(newcenterline);
    newVertexFilter->Update();
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection(newVertexFilter->GetOutputPort());
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper(newMapper);
    */

    // select the vertexes with z value higher than 0
    /*
    vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    for(vtkIdType i = 0; i<m_colon->GetOutput()->GetNumberOfPoints(); i++)
    {
        double p[3];
        m_colon->GetOutput()->GetPoint(i, p);
        //std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
        if(p[2] >= -213.051)  // the z value of the deformed centerline
        {
            ids->InsertNextValue(i);
        }
    }

    vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::POINT);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);
    selectionNode->GetProperties()->Set(vtkSelectionNode::CONTAINING_CELLS(), 1);

    vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
    extractSelection->SetInputData(0, m_colon->GetOutput());
    extractSelection->SetInputData(1, selection);
    extractSelection->Update();

    vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
    selected->ShallowCopy(extractSelection->GetOutput());
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(selected);
    geometryFilter->Update();
    */
    vtkSmartPointer<vtkPlane> clipPlane = vtkSmartPointer<vtkPlane>::New();
    clipPlane->SetOrigin(0, 0, -213.051);
    clipPlane->SetNormal(0, 0, 1);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(clipPlane);
    clipper->SetInputData(m_colon->GetOutput());

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputConnection(clipper->GetOutputPort());
    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    vtkSmartPointer<vtkFeatureEdges> featureEdges =
            vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->SetInputConnection(clipper->GetOutputPort());
    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> edgeMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    edgeMapper->SetInputConnection(featureEdges->GetOutputPort());
    vtkSmartPointer<vtkActor> edgeActor =
            vtkSmartPointer<vtkActor>::New();
    edgeActor->SetMapper(edgeMapper);
    edgeActor->GetProperty()->SetColor(1,0,0);

    m_showselectedwindow.show();
    m_showselectedwindow.RenderSelected(selectedActor);
    m_showselectedwindow.RenderSelected(edgeActor);


    /*
    double position[3] = {0,0, -30};
    vtkMath::Add(position, edgeActor->GetCenter(), position);
    m_showselectedwindow.GetRenderManager().GetLight()->SetPosition(position);

    double direction[3] = {0, 0, -1};
    vtkMath::Add(direction, edgeActor->GetCenter(), direction);
    m_showselectedwindow.GetRenderManager().GetLight()->SetFocalPoint(direction);

    vtkSmartPointer<vtkLightActor> lightActor = vtkSmartPointer<vtkLightActor>::New();
    lightActor->SetLight(m_showselectedwindow.GetRenderManager().GetLight());
    m_showselectedwindow.GetRenderManager().GetRender()->AddViewProp(lightActor);
    */





    //m_showselectedwindow.GetRenderManager().GetRender()->LightFollowCameraOff();

    m_showselectedwindow.GetRenderManager().addlight();
}
