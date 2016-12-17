#ifndef COLON_H
#define COLON_H
//ru
#include "object.h"

class Colon : public Object
{
private: int opacity;
public:
    Colon();
    ~Colon(){}
    void SetPoint(vtkSmartPointer<vtkPoints> points);
};

#endif // COLON_H
