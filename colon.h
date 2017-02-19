#ifndef COLON_H
#define COLON_H

#include "object.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>

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
};

#endif // COLON_H
