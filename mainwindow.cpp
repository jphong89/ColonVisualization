#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <unistd.h>

namespace
{
void CallbackFunction (vtkObject* caller, long unsigned int eventId,
                        void* clientData, void* callData );
void CallbackFunction_inverse (vtkObject* caller, long unsigned int eventId,
                        void* clientData, void* callData );
}
namespace
{
void CallbackFunction (vtkObject* caller,
                       long unsigned int vtkNotUsed(eventId),
                       void* clientData,
                       void* vtkNotUsed(callData) )
{
  vtkImageTracerWidget* tracerWidget = static_cast<vtkImageTracerWidget*>(caller);
  MainWindow* t_window = static_cast<MainWindow*>(clientData);

  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();

  tracerWidget->GetPath(path);

  if(path->GetNumberOfPoints() > 0)
  {
      double p[3];
      path->GetPoint(0, p);
      vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
      pointLocator->SetDataSet(t_window->GetData(0)->GetOutput());
      pointLocator->BuildLocator();
      int id = pointLocator->FindClosestPoint(p);
      std::cout<<id<<" :"<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;

      double p1[3], p2[3];
      t_window->GetData(0)->GetOutput()->GetPoint(id, p1);
      t_window->GetData(1)->GetOutput()->GetPoint(id, p2);
      vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
      points1->DeepCopy(t_window->GetTracerMark(0)->GetOutput()->GetPoints());
      points2->DeepCopy(t_window->GetTracerMark(1)->GetOutput()->GetPoints());
      points1->InsertNextPoint(p1);
      points2->InsertNextPoint(p2);

      t_window->GetTracerMark(0)->SetColor(0,0,1);
      t_window->GetTracerMark(1)->SetColor(0,0,1);
      t_window->GetTracerMark(0)->InputPoints(points1);
      t_window->GetTracerMark(1)->InputPoints(points2);
      t_window->GetRenderManager(0)->GetRender()->AddActor(t_window->GetTracerMark(0)->GetActor());
      //t_window->GetRenderManager(0)->renderModel(t_window->GetTracerMark(0)->GetActor());
      t_window->GetRenderManager(1)->GetRender()->AddActor(t_window->GetTracerMark(1)->GetActor());
      //t_window->GetRenderManager(1)->renderModel(t_window->GetTracerMark(1)->GetActor());

      QVTKWidget* left = t_window->findChild<QVTKWidget*>("left");
      left->GetRenderWindow()->Render();
      QVTKWidget* right = t_window->findChild<QVTKWidget*>("right");
      right->GetRenderWindow()->Render();
  }
}
void CallbackFunction_inverse(vtkObject* caller,
                       long unsigned int vtkNotUsed(eventId),
                       void* clientData,
                       void* vtkNotUsed(callData) )
{
  vtkImageTracerWidget* tracerWidget = static_cast<vtkImageTracerWidget*>(caller);
  MainWindow* t_window = static_cast<MainWindow*>(clientData);

  vtkSmartPointer<vtkPolyData> path =
    vtkSmartPointer<vtkPolyData>::New();

  tracerWidget->GetPath(path);

  if(path->GetNumberOfPoints() > 0)
  {
      double p[3];
      path->GetPoint(0, p);
      vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
      pointLocator->SetDataSet(t_window->GetData(1)->GetOutput());
      pointLocator->BuildLocator();
      int id = pointLocator->FindClosestPoint(p);

      double p1[3], p2[3];
      t_window->GetData(0)->GetOutput()->GetPoint(id, p1);
      t_window->GetData(1)->GetOutput()->GetPoint(id, p2);
      vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
      points1->DeepCopy(t_window->GetTracerMark(2)->GetOutput()->GetPoints());
      points2->DeepCopy(t_window->GetTracerMark(3)->GetOutput()->GetPoints());
      points1->InsertNextPoint(p1);
      points2->InsertNextPoint(p2);

      t_window->GetTracerMark(2)->SetColor(0,1,0);
      t_window->GetTracerMark(3)->SetColor(0,1,0);
      t_window->GetTracerMark(2)->InputPoints(points1);
      t_window->GetTracerMark(3)->InputPoints(points2);
      t_window->GetRenderManager(0)->GetRender()->AddActor(t_window->GetTracerMark(2)->GetActor());
      //t_window->GetRenderManager(0)->renderModel(t_window->GetTracerMark(2)->GetActor());
      t_window->GetRenderManager(1)->GetRender()->AddActor(t_window->GetTracerMark(3)->GetActor());
      //t_window->GetRenderManager(1)->renderModel(t_window->GetTracerMark(3)->GetActor());

      QVTKWidget* left = t_window->findChild<QVTKWidget*>("left");
      left->GetRenderWindow()->Render();
      QVTKWidget* right = t_window->findChild<QVTKWidget*>("right");
      right->GetRenderWindow()->Render();
  }
}
} // end anonymous namespace

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_filemanager = new FileManager;
    m_centerline = new Centerline;
    m_colon = new Colon;
    m_colon_new = new Colon;
    m_rendermanager = new RenderManager;
    m_rendermanager_right = new RenderManager;
    tracer = vtkSmartPointer<vtkImageTracerWidget>::New();
    tracer_inverse = vtkSmartPointer<vtkImageTracerWidget>::New();
    m_tracermark_1 = new TracerMarker;
    m_tracermark_2 = new TracerMarker;
    m_tracermark_inverse_1 = new TracerMarker;
    m_tracermark_inverse_2 = new TracerMarker;

    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->AddRenderer(m_rendermanager->GetRender());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->AddRenderer(m_rendermanager_right->GetRender());
}

MainWindow::~MainWindow()
{
    delete ui;
    delete m_filemanager;
    delete m_centerline;
    delete m_colon;
    delete m_colon_new;
    delete m_rendermanager;
    delete m_rendermanager_right;
    delete m_tracermark_1;
    delete m_tracermark_2;
    delete m_tracermark_inverse_1;
    delete m_tracermark_inverse_2;
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

    m_rendermanager->renderModel(m_colon->GetActor());
    QVTKWidget* left = this->findChild<QVTKWidget*>("left");

    tracer->GetLineProperty()->SetLineWidth(5);
    tracer->SetInteractor(left->GetRenderWindow()->GetInteractor());
    tracer->SetViewProp(m_colon->GetActor());
    tracer->AutoCloseOn();

    left->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction);
    callback->SetClientData(this);
    tracer->AddObserver(vtkCommand::EndInteractionEvent, callback);

    QVTKWidget* right = this->findChild<QVTKWidget*>("right");
    right->GetRenderWindow()->Render();
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

    // get the origin point
    /*
    double p[3];
    m_centerline->GetOutput()->GetPoint(0, p);
    std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<endl;
    */

    vtkSmartPointer<vtkPolyData> newColonPoly = vtkSmartPointer<vtkPolyData>::New();
    newColonPoly = m_centerline->EliminateTorsion(m_rendermanager,m_rendermanager_right, m_colon->GetOutput(), m_filemanager);
    m_colon_new->Object::SetInput(newColonPoly);

    m_filemanager->SaveFile(m_centerline->GetOutput(), "ModifiedCenterline.vtp");

    m_rendermanager->renderModel(m_centerline->GetActor());
    m_rendermanager_right->renderModel(m_colon_new->GetActor());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");

    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(m_colon_new->GetActor());

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();
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

void MainWindow::on_action_Deform_Colon_triggered(bool test)
{
    test = true;
    double factor = 2, r0 = 18.5793, adjust = 0.75;
    double k = 0.5, b;
    b = r0*(1-k);
    double aver = 0;
    double origin[3] = {667.6, 491.462, -213.051}, normal[3] = {0,0,1}, d1[3] = {1,0,0}, d2[3] = {0,1,0};

    if(test)
    {
        // configuration
        r0 = 2.79012; adjust = 0.9; b = r0*(1-k);
        double f1[3] = {7.904, -4.674, 11.167}, f2[3] = {-5.976, 3.461, -3.272}, f3[3] = {-0.826, -1.501, 10.889};
        origin[0] = f2[0];
        origin[1] = f2[1];
        origin[2] = f2[2];
        /* // visualize the fiducial points
        vtkSmartPointer<vtkPoints> fiducials = vtkSmartPointer<vtkPoints>::New();
        fiducials->InsertNextPoint(f1);
        fiducials->InsertNextPoint(f2);
        fiducials->InsertNextPoint(f3);
        vtkSmartPointer<vtkPolyData> fiducialspoly = vtkSmartPointer<vtkPolyData>::New();
        fiducialspoly->SetPoints(fiducials);
        vtkSmartPointer<vtkVertexGlyphFilter> vertexFiducial = vtkSmartPointer<vtkVertexGlyphFilter>::New();
        vertexFiducial->SetInputData(fiducialspoly);
        vertexFiducial->Update();
        vtkSmartPointer<vtkPolyDataMapper> fiducialsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        fiducialsMapper->SetInputConnection(vertexFiducial->GetOutputPort());
        fiducialsMapper->Update();
        vtkSmartPointer<vtkActor> fiducialActor = vtkSmartPointer<vtkActor>::New();
        fiducialActor->SetMapper(fiducialsMapper);
        fiducialActor->GetProperty()->SetPointSize(5);
        m_rendermanager->renderModel(fiducialActor);
        */

        double pd1[3], pd2[3], pd3[3], vpar[3], v2[3];
        double projection;
        vtkMath::Subtract(f1, f2, pd1);
        vtkMath::Normalize(pd1);
        vtkMath::Subtract(f3, f2, v2);
        projection = vtkMath::Dot(pd1, v2);
        vpar[0] = projection * pd1[0];
        vpar[1] = projection * pd1[1];
        vpar[2] = projection * pd1[2];
        vtkMath::Subtract(v2, vpar, pd2);

        //vtkMath::MultiplyScalar(pd2, -1);

        vtkMath::Normalize(pd2);
        vtkMath::Cross(pd1, pd2, pd3);
        vtkMath::Normalize(pd3);

        vtkSmartPointer<vtkPoints> newpoints = vtkSmartPointer<vtkPoints>::New();
        for(int i = 0; i < m_colon->GetOutput()->GetNumberOfPoints(); i++)
        {
            double p[3], v[3];
            m_colon->GetOutput()->GetPoint(i, p);
            vtkMath::Subtract(p, f2, v);
            double x, y, z;
            x = vtkMath::Dot(v, pd1);
            y = vtkMath::Dot(v, pd2);
            z = vtkMath::Dot(v, pd3);
            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, y);
            vtkMath::MultiplyScalar(vz, z);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            newpoints->InsertNextPoint(newp);
        }
        m_colon->GetOutput()->SetPoints(newpoints);
        m_rendermanager->renderModel(m_colon->GetActor());
    }

    vtkSmartPointer<vtkPlane> clipPlane = vtkSmartPointer<vtkPlane>::New();
    clipPlane->SetOrigin(origin);
    clipPlane->SetNormal(normal);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(clipPlane);
    clipper->SetInputData(m_colon->GetOutput());
    clipper->Update();
    vtkSmartPointer<vtkPolyData> leftpart = vtkSmartPointer<vtkPolyData>::New();
    leftpart->DeepCopy(clipper->GetOutput());
    vtkSmartPointer<vtkPoints> leftnewpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i < leftpart->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        leftpart->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        double x = vtkMath::Dot(d1, v);
        double z = vtkMath::Dot(normal, v);
        if(z < 0) z = -z;
        double y = vtkMath::Dot(d2, v);
        double r = sqrt(z*z + y*y);
        double angle = acos(y/r);
        //double newr = r * factor * exp(adjust *((r0 - r)/r0));
        double newr = factor * (r * k + b);
        double newangle = angle / factor;

        double newy = newr * cos(newangle);
        double newz = newr * sin(newangle);

        double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
        vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
        vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
        vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
        vtkMath::MultiplyScalar(vx, x);
        vtkMath::MultiplyScalar(vy, newy);
        vtkMath::MultiplyScalar(vz, newz);
        vtkMath::Add(origin, vx, temp1);
        vtkMath::Add(temp1, vy, temp2);
        vtkMath::Add(temp2, vz, newp);
        leftnewpoints->InsertNextPoint(newp);
        aver += r;
            //std::cout<<newp[0]<<" "<<newp[1]<<" "<<newp[2]<<endl;
    }
    leftpart->SetPoints(leftnewpoints);

    vtkSmartPointer<vtkPolyData> rightpart = vtkSmartPointer<vtkPolyData>::New();
    vtkMath::MultiplyScalar(normal, -1);
    clipPlane->SetNormal(normal);
    clipper->SetClipFunction(clipPlane);
    clipper->Update();
    rightpart->DeepCopy(clipper->GetOutput());
    vtkSmartPointer<vtkPoints> rightnewpoints = vtkSmartPointer<vtkPoints>::New();
    for(int i=0; i < rightpart->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        rightpart->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        double x = vtkMath::Dot(d1, v);
        double z = vtkMath::Dot(normal, v);
        if(z < 0) z = -z;
        double y = vtkMath::Dot(d2, v);
        double r = sqrt(z*z + y*y);
        double angle = acos(y/r);
        //double newr = r * factor * exp(adjust *((r0 - r)/r0));
        double newr = factor * (r * k + b);
        double newangle = angle / factor;

        double newy = newr * cos(newangle);
        double newz = newr * sin(newangle);

        double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
        vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
        vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
        vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
        vtkMath::MultiplyScalar(vx, x);
        vtkMath::MultiplyScalar(vy, newy);
        vtkMath::MultiplyScalar(vz, newz);
        vtkMath::Add(origin, vx, temp1);
        vtkMath::Add(temp1, vy, temp2);
        vtkMath::Add(temp2, vz, newp);
        rightnewpoints->InsertNextPoint(newp);
        aver += r;
            //std::cout<<newp[0]<<" "<<newp[1]<<" "<<newp[2]<<endl;
    }
    rightpart->SetPoints(rightnewpoints);
    aver = aver / (leftpart->GetNumberOfPoints() + rightpart->GetNumberOfPoints());
    std::cout<<"Average radius is "<<aver<<endl;

    // one to one mapping used for interaction
    vtkSmartPointer<vtkPoints> onetoonepoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> onetoone = vtkSmartPointer<vtkPolyData>::New();
    onetoone->DeepCopy(m_colon->GetOutput());
    for(vtkIdType i=0; i<m_colon->GetOutput()->GetNumberOfPoints(); i++)
    {
        double p[3], v[3];
        m_colon->GetOutput()->GetPoint(i, p);
        vtkMath::Subtract(p, origin, v);
        if(vtkMath::Dot(v, normal) > 0) // left points
        {
            double x = vtkMath::Dot(d1, v);
            double z = vtkMath::Dot(normal, v);
            if(z < 0) z = -z;
            double y = vtkMath::Dot(d2, v);
            double r = sqrt(z*z + y*y);
            double angle = acos(y/r);
            //double newr = r * factor * exp(adjust *((r0 - r)/r0));
            double newr = factor * (r * k + b);
            double newangle = angle / factor;

            double newy = newr * cos(newangle);
            double newz = newr * sin(newangle);

            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, newy);
            vtkMath::MultiplyScalar(vz, newz);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            onetoonepoints->InsertNextPoint(newp);
        }
        else // right points
        {
            vtkMath::MultiplyScalar(normal, -1);
            double x = vtkMath::Dot(d1, v);
            double z = vtkMath::Dot(normal, v);
            if(z < 0) z = -z;
            double y = vtkMath::Dot(d2, v);
            double r = sqrt(z*z + y*y);
            double angle = acos(y/r);
            //double newr = r * factor * exp(adjust *((r0 - r)/r0));
            double newr = factor * (r * k + b);
            double newangle = angle / factor;

            double newy = newr * cos(newangle);
            double newz = newr * sin(newangle);

            double vx[3], vy[3], vz[3], temp1[3], temp2[3], newp[3];
            vx[0] = d1[0]; vx[1] = d1[1]; vx[2] = d1[2];
            vy[0] = d2[0]; vy[1] = d2[1]; vy[2] = d2[2];
            vz[0] = normal[0]; vz[1] = normal[1]; vz[2] = normal[2];
            vtkMath::MultiplyScalar(vx, x);
            vtkMath::MultiplyScalar(vy, newy);
            vtkMath::MultiplyScalar(vz, newz);
            vtkMath::Add(origin, vx, temp1);
            vtkMath::Add(temp1, vy, temp2);
            vtkMath::Add(temp2, vz, newp);
            onetoonepoints->InsertNextPoint(newp);
            vtkMath::MultiplyScalar(normal, -1);
        }
    }
    onetoone->SetPoints(onetoonepoints);

    vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputData(leftpart);
    appendFilter->AddInputData(rightpart);
    appendFilter->Update();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    vtkSmartPointer<vtkPolyData> bended = vtkSmartPointer<vtkPolyData>::New();
    bended->DeepCopy(cleanFilter->GetOutput());

    vtkSmartPointer<vtkPolyDataMapper> selectedMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    selectedMapper->SetInputData(bended);
    //selectedMapper->SetInputData(onetoone);
    vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
    selectedActor->SetMapper(selectedMapper);

    vtkSmartPointer<vtkFeatureEdges> featureEdges =
            vtkSmartPointer<vtkFeatureEdges>::New();
    featureEdges->SetInputData(bended);
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
    edgeActor->GetProperty()->SetColor(255,0,0);

    m_colon_new->Object::SetInput(onetoone);
    m_rendermanager_right->renderModel(selectedActor);
    //m_rendermanager_right->renderModel(edgeActor);
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");

    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(selectedActor);

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();

    //m_showselectedwindow.show();
    //m_showselectedwindow.RenderSelected(selectedActor);

    m_filemanager->SaveFile(bended, "openedcolontextured.off", true);

    //m_showselectedwindow.RenderSelected(edgeActor);
    //m_showselectedwindow.GetRenderManager().GetRender()->SetBackground(0.1, 0.6, 1);
}

void MainWindow::on_tracer_toggled(bool checked)
{
    if(checked)
    {
        tracer->On();
    }
    else
    {
        tracer->Off();

        m_rendermanager->GetRender()->RemoveActor(m_tracermark_1->GetActor());
        m_rendermanager_right->GetRender()->RemoveActor(m_tracermark_2->GetActor());
        vtkSmartPointer<vtkPoints> empty = vtkSmartPointer<vtkPoints>::New();
        m_tracermark_1->InputPoints(empty);
        m_tracermark_2->InputPoints(empty);

        QVTKWidget* left = this->findChild<QVTKWidget*>("left");
        left->GetRenderWindow()->Render();
        QVTKWidget* right = this->findChild<QVTKWidget*>("right");
        right->GetRenderWindow()->Render();
    }
}

void MainWindow::on_action_Load_Deformed_Surface_triggered()
{
    //Open file dialog
    QString filePath = QFileDialog::getOpenFileName(
                this, tr("Open File"), "",
                tr("3Dmodels (*.stl *.vtp *.ply *.obj *.off)"));
    if(filePath.isEmpty()) return;

    m_filemanager->LoadNewFile(filePath);
    m_colon_new->Object::SetInput(m_filemanager->getfile());

    m_rendermanager_right->renderModel(m_colon_new->GetActor());
    QVTKWidget* right = this->findChild<QVTKWidget*>("right");

    tracer_inverse->GetLineProperty()->SetLineWidth(5);
    tracer_inverse->SetInteractor(right->GetRenderWindow()->GetInteractor());
    tracer_inverse->SetViewProp(m_colon_new->GetActor());

    right->GetRenderWindow()->Render();
    vtkSmartPointer<vtkCallbackCommand> callback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    callback->SetCallback(CallbackFunction_inverse);
    callback->SetClientData(this);
    tracer_inverse->AddObserver(vtkCommand::EndInteractionEvent, callback);


    QVTKWidget* left = this->findChild<QVTKWidget*>("left");
    left->GetRenderWindow()->Render();
}

void MainWindow::on_checkBox_toggled(bool checked)
{
    if(checked)
    {
        tracer_inverse->On();
    }
    else
    {
        tracer_inverse->Off();

        m_rendermanager->GetRender()->RemoveActor(m_tracermark_inverse_1->GetActor());
        m_rendermanager_right->GetRender()->RemoveActor(m_tracermark_inverse_2->GetActor());
        vtkSmartPointer<vtkPoints> empty = vtkSmartPointer<vtkPoints>::New();
        m_tracermark_inverse_1->InputPoints(empty);
        m_tracermark_inverse_2->InputPoints(empty);

        QVTKWidget* left = this->findChild<QVTKWidget*>("left");
        left->GetRenderWindow()->Render();
        QVTKWidget* right = this->findChild<QVTKWidget*>("right");
        right->GetRenderWindow()->Render();
    }
}
