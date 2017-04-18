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

#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkImageActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkImageTracerWidget.h>
#include <vtkImageMapper3D.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkInteractorStyleImage.h>
#include <vtkProperty.h>

#include "tracermarker.h"

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
    RenderManager* GetRenderManager(int side = 0)
    {
        if(side == 0) return m_rendermanager;
        else return m_rendermanager_right;
    }
    Colon* GetData(int side = 0)
    {
        if(side == 0) return m_colon;
        else return m_colon_new;
    }
    TracerMarker* GetTracerMark(int id = 0)
    {
        if(id == 0)
            return m_tracermark_1;
        else if(id == 1)
            return m_tracermark_2;
        else if(id == 2)
            return m_tracermark_inverse_1;
        else
            return m_tracermark_inverse_2;
    }

private slots:
    void on_actionNew_file_triggered();

    void on_actionLoad_Centerline_triggered();

    void on_action_Deform_Colon_triggered(bool test = false);

    void on_tracer_toggled(bool checked);

    void on_action_Load_Deformed_Surface_triggered();

    void on_checkBox_toggled(bool checked);

private:
    Ui::MainWindow *ui;
    FileManager * m_filemanager;
    Centerline *m_centerline;
    Colon *m_colon;
    Colon *m_colon_new;
    RenderManager *m_rendermanager;
    RenderManager *m_rendermanager_right;
    QTimer *m_timer;

    vtkSmartPointer<vtkImageTracerWidget> tracer;
    vtkSmartPointer<vtkImageTracerWidget> tracer_inverse;
    TracerMarker *m_tracermark_1;
    TracerMarker *m_tracermark_2;
    TracerMarker *m_tracermark_inverse_1;
    TracerMarker *m_tracermark_inverse_2;

    ShowSelectedWindow m_showselectedwindow;
};

#endif // MAINWINDOW_H
