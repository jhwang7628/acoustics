#include <ui/FDTD_AcousticSimulator_Widget.h>

//##############################################################################
//##############################################################################
FDTD_AcousticSimulator_Widget::
FDTD_AcousticSimulator_Widget(std::shared_ptr<FDTD_AcousticSimulator_Viewer> &viewer, const RunMode &runMode)
    : _viewer(viewer), _solverSettings(_viewer->_solverSettings), _runMode(runMode)
{
    // initialize layout
    _layout                = new QGridLayout(); 
    _slider_impulseScaling = new QSlider(Qt::Horizontal, this); 
    _slider_impulseScaling->setRange(-3./IMP_SLIDER_SCALE, 3./IMP_SLIDER_SCALE); 
    _slider_impulseScaling->setValue(0);
     _label_impulseScaling = new QLabel("0"); 
     _label_impulseScaling->setFixedWidth(100);
      _text_impulseScaling = new QLabel("Impulse Scaling"); 
      _text_impulseScaling->setFixedWidth(100);
    _slider_simulationTimeline = new QSlider(Qt::Horizontal, this); 
    _slider_simulationTimeline->setRange(0, (int)ceil(_solverSettings->timeEnd/_solverSettings->timeStepSize)); 
    _slider_simulationTimeline->setValue(0);
     _label_simulationTimeline = new QLabel("0"); 
     _label_simulationTimeline->setFixedWidth(200);
      _text_simulationTimeline = new QLabel("Simulation Time"); 
      _text_simulationTimeline->setFixedWidth(200);
    _layout->addWidget(_viewer.get()         , 0, 0, 1, 3);
    _layout->addWidget(  _text_impulseScaling, 1, 0);
    _layout->addWidget(_slider_impulseScaling, 1, 1);
    _layout->addWidget( _label_impulseScaling, 1, 2);
    _layout->addWidget(  _text_simulationTimeline, 2, 0);
    _layout->addWidget(_slider_simulationTimeline, 2, 1);
    _layout->addWidget( _label_simulationTimeline, 2, 2);
    setLayout(_layout);
    resize(800, 600);
    // signal-slot stuff
    connect(_slider_impulseScaling, SIGNAL(valueChanged(int)), this, SLOT(SliderValueChanged()));
    connect(_slider_simulationTimeline, SIGNAL(valueChanged(int)), this, SLOT(SliderValueChanged()));
}

//##############################################################################
//##############################################################################
FDTD_AcousticSimulator_Widget::
~FDTD_AcousticSimulator_Widget()
{
    delete _layout; // this object should have ownerships to all its widget children and will delete them
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Widget::
SliderValueChanged()
{
    const REAL impulseScaling = (REAL)_slider_impulseScaling->value()*IMP_SLIDER_SCALE; 
    const REAL newTime = (REAL)_slider_simulationTimeline->value()*_solverSettings->timeStepSize; 
    _label_impulseScaling->setText(QString("%1").arg(impulseScaling));
    _label_simulationTimeline->setText(QString("%1").arg(newTime));
    // update viewer
    _viewer->_drawImpulseScaling = pow(10., impulseScaling); 
    _viewer->updateGL();
}
