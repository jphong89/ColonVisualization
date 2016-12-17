#include "lightdialog.h"
#include "genesyndata.h"
#include "ui_lightdialog.h"
#include <QVTKWidget.h>
#include "mainwindow.h"


lightDialog::lightDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::lightDialog)
{
    ui->setupUi(this);
}



void lightDialog::setscene(RenderManager *t_scene)
{
    m_scene = t_scene;
}

void lightDialog::setwindow(MainWindow *t_window)
{
    m_window = t_window;
}

lightDialog::~lightDialog()
{
    delete ui;
}

void lightDialog::on_radioButton_5_clicked()
{
    m_scene->setconstantlight(this->findChild<QSlider*>("Attenvalue")->value());
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

void lightDialog::on_radioButton_3_clicked()
{
    m_scene->setlinearlight(this->findChild<QSlider*>("Attenvalue")->value());
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

void lightDialog::on_radioButton_4_clicked()
{
    m_scene->setqudraticlight(this->findChild<QSlider*>("Attenvalue")->value());
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

void lightDialog::on_horizontalSlider_valueChanged(int value)
{
    m_scene->setambient(value);
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

void lightDialog::on_horizontalSlider_2_valueChanged(int value)
{
    m_scene->setdiffuse(value);
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}

void lightDialog::on_horizontalSlider_3_valueChanged(int value)
{
    m_scene->setspecular(value);
    QVTKWidget* widget = m_window->findChild<QVTKWidget*>("qvtk");
    widget->GetRenderWindow()->Render();
}
void lightDialog::on_checkBox_clicked(bool checked)
{
    if(checked) m_window->addlight();
    else m_scene->GetRender()->RemoveAllLights();
}
