#ifndef FILEMANAGER_H
#define FILEMANAGER_H

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <QString>
#include <vtkSTLReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>

class FileManager
{
public:
    FileManager();
    ~FileManager();
    void LoadNewFile(QString fileName);
    void SaveFile(vtkSmartPointer<vtkPolyData> polydata, char *filename);
    vtkSmartPointer<vtkPolyData> getfile() {return m_polydata;}

private:
    vtkSmartPointer<vtkXMLPolyDataReader> m_vtp;
    vtkSmartPointer<vtkSTLReader> m_stl;
    vtkSmartPointer<vtkOBJReader> m_obj;
    vtkSmartPointer<vtkPLYReader> m_ply;

    vtkSmartPointer<vtkPolyData> m_polydata;
};

#endif // FILEMANAGER_H