#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include "centerline.h"
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void constructA(int idx, std::map<int, double> &coefficientMap);
void constructb(int idx, double *b);

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed);

void updateA(int m, int i1, int i2, double weight,std::map<int, double> &coefficientMap);

#endif // OPTIMIZATION_H
