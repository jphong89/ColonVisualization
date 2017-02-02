#include "filemanager.h"

FileManager::FileManager()
{
    m_vtp = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    m_stl = vtkSmartPointer<vtkSTLReader>::New();
    m_obj = vtkSmartPointer<vtkOBJReader>::New();
    m_ply = vtkSmartPointer<vtkPLYReader>::New();

    m_polydata = vtkSmartPointer<vtkPolyData>::New();
}
FileManager::~FileManager()
{
}

void FileManager::LoadNewFile(QString fileName) // Read 4 different file formats
{
    if(fileName.endsWith("vtp"))
    {
        m_vtp->SetFileName(fileName.toStdString().c_str());
        m_vtp->Update();
        m_polydata = m_vtp->GetOutput();
    }
    else if(fileName.endsWith("stl"))
    {
        m_stl->SetFileName(fileName.toStdString().c_str());
        m_stl->Update();
        m_polydata = m_stl->GetOutput();
    }
    else if(fileName.endsWith("obj"))
    {
        m_obj->SetFileName(fileName.toStdString().c_str());
        m_obj->Update();
        m_polydata = m_obj->GetOutput();
    }
    else if(fileName.endsWith("ply"))
    {
        m_ply->SetFileName(fileName.toStdString().c_str());
        m_ply->Update();
        m_polydata = m_ply->GetOutput();
    }
    else if(fileName.endsWith("off"))
    {
        m_polydata->DeepCopy(readFileOFF(fileName.toStdString().c_str()));
    }
}

vtkSmartPointer<vtkPolyData> FileManager::readFileOFF(const char *filename)
{
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    FILE *in = fopen(filename, "r");
    if(!in)
    {
        std::cerr<<"failed to open file!"<<endl;
        exit(1);
    }
    int n_pts, n_polys;
    fscanf(in, "COFF\n%d %d 0\n", &n_pts, &n_polys);
    std::cout<<"points: "<<n_pts<<" polys: "<<n_polys<<endl;
    for(int i = 0; i<n_pts; i++)
    {
        float p[3];
        unsigned char c[3];
        fscanf(in, "%f %f %f %hhu %hhu %hhu 255\n", &p[0], &p[1], &p[2], &c[0], &c[1], &c[2]);
        points->InsertNextPoint(p);
        colors->InsertNextTypedTuple(c);
    }
    for(int i = 0; i<n_polys; i++)
    {
        vtkIdType vertices[3];
        fscanf(in, "3 %lld %lld %lld\n", &vertices[0], &vertices[1], &vertices[2]);
        //std::cout<<vertices[0]<<" "<<vertices[1]<<" "<<vertices[2]<<endl;
        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
        triangle->GetPointIds()->SetId(0, vertices[0]);
        triangle->GetPointIds()->SetId(1, vertices[1]);
        triangle->GetPointIds()->SetId(2, vertices[2]);
        polys->InsertNextCell(triangle);
    }
    output->SetPoints(points);
    output->SetPolys(polys);
    output->GetPointData()->SetScalars(colors);
    std::cout<<output->GetNumberOfPoints()<<" "<<output->GetNumberOfPolys()<<" "<<endl;
    return output;
}

void FileManager::SaveFile(vtkSmartPointer<vtkPolyData> polydata, char* filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(polydata);
    writer->Write();
}
