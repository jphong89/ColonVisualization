#ifndef COLON_H
#define COLON_H

#include "object.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <vtkFillHolesFilter.h>
#include <vtkPolyDataConnectivityFilter.h>

class Colon : public Object
{
private: int opacity;
public:
    Colon();
    ~Colon(){}
    void SetPoint(vtkSmartPointer<vtkPoints> points);
    void AddTexture();
    void SmoothSurface();
    void Decimation();
    void testDeformation();
    void FillHoles(vtkSmartPointer<vtkPolyData> object = NULL);
    void RemoveUnconnectedBlobs(vtkSmartPointer<vtkPolyData> object = NULL);
};

#endif // COLON_H
