#ifndef FDTD_ACOUSTIC_SIMULATOR_WIDGET_H
#define FDTD_ACOUSTIC_SIMULATOR_WIDGET_H
#include <memory>
#include <string>
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <QGridLayout>
#include <QSlider>
#include <QLabel>
#include <QtGui>

//##############################################################################
// Constants
//##############################################################################
#define IMP_SLIDER_SCALE 0.1


//##############################################################################
// Class that manage GUI interface including viewer.
//##############################################################################
class FDTD_AcousticSimulator_Widget : public QWidget
{
    Q_OBJECT
    private: 
        std::shared_ptr<FDTD_AcousticSimulator_Viewer> _viewer;
        QGridLayout *_layout;
        QSlider     *_slider_impulseScaling;
        QLabel      * _label_impulseScaling;

    private slots: 
        void SliderValueChanged(); 

    public: 
        FDTD_AcousticSimulator_Widget(std::shared_ptr<FDTD_AcousticSimulator_Viewer> &viewer); 
};

#endif
