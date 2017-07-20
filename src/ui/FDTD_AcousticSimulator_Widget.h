#ifndef FDTD_ACOUSTIC_SIMULATOR_WIDGET_H
#define FDTD_ACOUSTIC_SIMULATOR_WIDGET_H
#include <memory>
#include <string>
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <wavesolver/PML_WaveSolver_Settings.h>
#include <ui/QtHelper.h>

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
    public: 
        enum RunMode{NORMAL=0, PREVIEW=1};
        FDTD_AcousticSimulator_Widget(std::shared_ptr<FDTD_AcousticSimulator_Viewer> &viewer, const RunMode &runMode=NORMAL); 
        ~FDTD_AcousticSimulator_Widget();

    private: 
        std::shared_ptr<FDTD_AcousticSimulator_Viewer>  _viewer;
        std::shared_ptr<PML_WaveSolver_Settings>        _solverSettings;
        RunMode      _runMode;
        QGridLayout *_layout;
        QGridLayout *_controlPanelLayout;
        QSlider     *_slider_impulseScaling;
        QLabel      * _label_impulseScaling;
        QLabel      *  _text_impulseScaling;
        QSlider     *_slider_simulationTimeline;
        QLabel      * _label_simulationTimeline;
        QLabel      *  _text_simulationTimeline;
        QPushButton *_button_resetSimulation;
        QPushButton *_button_generateSlice_x;
        QPushButton *_button_generateSlice_y;
        QPushButton *_button_generateSlice_z;
        QPushButton *_button_clearPressures; 
        QPushButton *_button_clearSources; 

    private slots: 
        void SliderValueChanged(); 
        void ResetSystemTime(const bool &fromSlider = true);
        void GenerateSlice(const int &dim, const REAL &offset);
        void ClearSources(); 
};

#endif
