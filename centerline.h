#ifndef CENTERLINE_H
#define CENTERLINE_H
#include "object.h"
#include "rendermanager.h"
#include "filemanager.h"
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <vtkPlane.h>
#include <vtkPlanes.h>
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
#include <vtkMatrix3x3.h>
#include <vtkTriangleFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>


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
    void SmoothCenterline(double std, vtkSmartPointer<vtkIntArray> ViolationNums = NULL);
    void Smooth(vtkSmartPointer<vtkDoubleArray> Candidate, double std);
    void GaussianTangents(vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void GaussianNormals(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures, vtkSmartPointer<vtkDoubleArray> Tangents, double std);
    void UniformSample(int resolution);
    void splineTangent(double *tangent, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);
    void splineNormal(double *normal, vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize, double &curvature);
    void PutNormalsOnSameSide(vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Curvatures);
    double splineTorsion(vtkSmartPointer<vtkParametricSpline> spline, double t_u, double stepsize);

    vtkSmartPointer<vtkPolyData> EliminateTorsion(RenderManager* t_rendermanager, vtkSmartPointer<vtkPolyData> t_colon);
    vtkSmartPointer<vtkPolyData> Deformation(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                                             vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals,
                                             vtkSmartPointer<vtkPolyData> t_colon, RenderManager *t_rendermanager,
                                             vtkSmartPointer<vtkDoubleArray> PlaneOriginals, vtkSmartPointer<vtkDoubleArray> PlaneNormals,
                                             vtkSmartPointer<vtkDoubleArray> RefDirections);
    void VisualizeTNB(vtkSmartPointer<vtkDoubleArray> S, vtkSmartPointer<vtkDoubleArray> Curvatures,
                      vtkSmartPointer<vtkDoubleArray> Tangents, vtkSmartPointer<vtkDoubleArray> Normals, vtkSmartPointer<vtkDoubleArray> Binormals,
                      RenderManager* t_rendermanager);
    void VisualizeSpoke(vtkSmartPointer<vtkPoints> CurvaturePoints, vtkSmartPointer<vtkIntArray> ViolationNums, RenderManager *t_rendermanager);
    void VisualizeOriginalCurve(RenderManager* t_rendermanager);

    vtkSmartPointer<vtkPolyData> ReorderContour(vtkSmartPointer<vtkPolyData> cutCircle);
    double SinglePath(double **costVrt, double **costHrz,int sx, int sy, int tx, int ty, int *steps);
    void ContoursToSurface(RenderManager* t_rendermanager, FileManager* t_filemanager);
};

#endif // CENTERLINE_H
