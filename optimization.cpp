#include "optimization.h"

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed)
{
    vtkSmartPointer<vtkPolyData> OptimizedSurface = vtkSmartPointer<vtkPolyData>::New();
    OptimizedSurface->DeepCopy(SurfaceLineUp);
    // maintain a list to record the correspondence between original ids and the variable ids
    vtkSmartPointer<vtkIdList> Ids = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(!Is_Fixed[i])
        {
            Ids->InsertNextId(i);
        }
    }
    int m = Ids->GetNumberOfIds() * 3;
    std::cout<<"total number of variables = "<<m<<endl;
    std::map<int, double> coefficientMap;
    coefficientMap.clear();
    return OptimizedSurface;
}

void updateA(int m, int i1, int i2, double weight,std::map<int, double> &coefficientMap)
{
    std::map<int, double>::iterator iter;
    int k = (i1*3) + (i2*3) * m;

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<int,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*3+1) + (i2*3+1) * m;

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<int,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*3+2) + (i2*3+2) * m;

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<int,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }
}

void constructA(int idx, std::map<int, double> &coefficientMap)
{

}
void constructb(int idx, double *b)
{

}
