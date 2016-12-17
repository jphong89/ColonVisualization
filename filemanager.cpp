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
}

void FileManager::SaveFile(vtkSmartPointer<vtkPolyData> polydata, char* filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(polydata);
    writer->Write();
}
