#include "colon.h"
#include <iostream>
Colon::Colon()
{
    actor->GetProperty()->SetOpacity(1);
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
void Colon::SmoothSurface()
{
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(model);
    smoothFilter->SetNumberOfIterations(15);
    smoothFilter->SetRelaxationFactor(0.1);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    model->DeepCopy(smoothFilter->GetOutput());
}
void Colon::Decimation()
{
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
      decimate->SetInputData(model);
      decimate->SetTargetReduction(.8); //10% reduction (if there was 100 triangles, now there will be 90)
      decimate->Update();
      model->DeepCopy(decimate->GetOutput());
}
void Colon::testDeformation()
{
    ofstream file;
    file.open("/home/ruibinma/Desktop/testDeformation.txt");
    file<<100<<" "<<100<<" "<<100<<" "<<1<<std::endl;
    for(vtkIdType i = 1; i < model->GetNumberOfPoints(); i++)
    {
        file<<1<<" "<<1<<" "<<1<<" "<<0<<std::endl;
    }
    file.close();
}
