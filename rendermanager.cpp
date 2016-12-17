#include "rendermanager.h"

RenderManager::RenderManager()
{
    render = vtkSmartPointer<vtkRenderer>::New();

    light = vtkSmartPointer<vtkLight>::New();
    light->SetLightTypeToHeadlight();
    light->SetAttenuationValues(0, 0.01, 0);
    light->SetPositional(true);
    light->SetConeAngle(180);
    light->SetDiffuseColor(10,10,10);
    light->SetAmbientColor(10,10,10);
    light->SetSpecularColor(10,10,10);
    light->SetLightTypeToHeadlight();
}

void RenderManager::renderModel(vtkSmartPointer<vtkActor> t_actor)
{
    render->AddActor(t_actor);
    render->ResetCamera();
}

void RenderManager::addlight()
{
    render->AddLight(light);
}
void RenderManager::setconstantlight(int value)
{
    light->SetAttenuationValues(float(value)/100,0,0);
}

void RenderManager::setlinearlight(int value)
{
    light->SetAttenuationValues(0,float(value)/100,0);
}

void RenderManager::setqudraticlight(int value)
{
    light->SetAttenuationValues(0,0,float(value)/100);
}

void RenderManager::setambient(int value)
{
    light->SetAmbientColor(value,value,value);
}

void RenderManager::setdiffuse(int value)
{
    light->SetDiffuseColor(value,value,value);
}

void RenderManager::setspecular(int value)
{
    light->SetSpecularColor(value,value,value);
}
