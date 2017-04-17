#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFileDialog>
#include "filemanager.h"
#include "centerline.h"
#include "object.h"
#include "colon.h"
#include "rendermanager.h"
#include "showselectedwindow.h"
#include "genesyndata.h"
#include <QVTKWidget.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkCutter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkDijkstraImageGeodesicPath.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkInformation.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkGeometryFilter.h>
#include <QTimer>
#include <vtkLightCollection.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkCardinalSpline.h>
#include <vtkSplineFilter.h>
#include <vtkClipPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkLightActor.h>
#include "vtkpolydatagroup.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void addlight();
    vtkSmartPointer<vtkPolyData> GeodesicPath(vtkSmartPointer<vtkPolyData> points);
    vtkSmartPointer<vtkPolyData> Upsampling(vtkSmartPointer<vtkPolyData> path);

private slots:
    void on_actionNew_file_triggered();

    void on_actionLoad_Centerline_triggered();

    void on_action_Deform_Colon_triggered(bool test = false);

private:
    Ui::MainWindow *ui;
    FileManager * m_filemanager;
    Centerline *m_centerline;
    Colon *m_colon;
    RenderManager *m_rendermanager;
    QTimer *m_timer;

    ShowSelectedWindow m_showselectedwindow;
};

#endif // MAINWINDOW_H
