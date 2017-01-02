#include "colon.h"
#include <iostream>
Colon::Colon()
{
    actor->GetProperty()->SetOpacity(0.5);
    //std::cout<<"Colon Opacity: "<<actor->GetProperty()->GetOpacity()<<endl;
}
void Colon::SetPoint(vtkSmartPointer<vtkPoints> points)
{
    model->SetPoints(points);
    mapper->Update();
}
