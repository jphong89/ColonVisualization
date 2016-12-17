#ifndef CENTERLINE_H
#define CENTERLINE_H
#include "object.h"
#include "rendermanager.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <vtkPlane.h>
#include <vtkPlaneSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSplineFilter.h>
#include <vtkCardinalSpline.h>
#include <vtkGlyph3D.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkThinPlateSplineTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPointLocator.h>
#include <vtkCutter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>

class Centerline : public Object
{
private:
    double* Gaussian1DKernel(int len, double std);
    double* GaussianDerivative1DKernel(int len, double std);
    double* conv(double *A, double *B, int lenA, int lenB, int *lenC);
public:
    Centerline();
    ~Centerline(){}
    int GetNumberOfPoints();
    vtkSmartPointer<vtkPlane> GetVerticalPlane(vtkIdType i, double* normalDirection);
    void SmoothCenterline(double std);
    void GaussianTangents(vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void GaussianNormals(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void UniformSample(int resolution);
    void splineTangent(double *tangent, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);
    void splineNormal(double *normal, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize, double &curvature);
    void PutNormalsOnSameSide(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures);
    double splineTorsion(vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);

    vtkSmartPointer<vtkPolyData> EliminateTorsion(RenderManager* t_rendermanager, vtkSmartPointer<vtkPolyData> t_colon);
    vtkSmartPointer<vtkPolyData> Deformation(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                                             RenderManager* t_rendermanager, vtkSmartPointer<vtkPolyData> t_colon);
    void VisualizeTNB(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                      vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                      RenderManager* t_rendermanager);
};

#endif // CENTERLINE_H