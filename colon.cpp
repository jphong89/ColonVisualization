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
void Colon::AddTexture()
{
    unsigned char color[3] = {255, 200, 200};
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetName("Color");
    colors->SetNumberOfComponents(3);
    for(vtkIdType i=0; i<model->GetNumberOfPoints(); i++)
    {
        colors->InsertNextTypedTuple(color);
    }
    model->GetPointData()->SetScalars(colors);
}
