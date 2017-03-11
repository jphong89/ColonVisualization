#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "centerline.h"
#include <Eigen/Sparse>
#include <vtkExtractEdges.h>

#define REGWEIGHT 1
#define BENDWEIGHT 1

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void constructA(std::map<int, double> &coefficientMap,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp,
                bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids,int* InvertIds);
void constructb(int idx, double *b,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids);

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed);

void updateA(int m, int i1, int i2, double weight,std::map<int, double> &coefficientMap);

vtkSmartPointer<vtkIdList> GetFacetsOfEdge(vtkSmartPointer<vtkPolyData> mesh, int idx1, int idx2);

double ComputeStretchWeight(vtkSmartPointer<vtkPolyData> t_colon, int idx1, int idx2);

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id);

double ComputeFacetArea(vtkSmartPointer<vtkPolyData> t_colon,  vtkIdType cellid);
#endif // OPTIMIZATION_H
