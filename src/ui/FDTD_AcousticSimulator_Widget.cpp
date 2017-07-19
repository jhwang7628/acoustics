#include <ui/FDTD_AcousticSimulator_Widget.h>

//##############################################################################
//##############################################################################
FDTD_AcousticSimulator_Widget::
FDTD_AcousticSimulator_Widget(std::shared_ptr<FDTD_AcousticSimulator_Viewer> &viewer, const RunMode &runMode)
    : _viewer(viewer), _solverSettings(_viewer->_solverSettings), _runMode(runMode)
{
    // initialize layout
    _layout                = new QGridLayout(); 
    _controlPanelLayout    = new QGridLayout(); 
    _slider_impulseScaling = new QSlider(Qt::Horizontal, this); 
    _slider_impulseScaling->setRange(-3./IMP_SLIDER_SCALE, 3./IMP_SLIDER_SCALE); 
    _slider_impulseScaling->setValue(0);
     _label_impulseScaling = new QLabel("0"); 
     _label_impulseScaling->setFixedWidth(100);
      _text_impulseScaling = new QLabel("Impulse Scaling"); 
      _text_impulseScaling->setFixedWidth(100);
    _slider_simulationTimeline = new QSlider(Qt::Horizontal, this); 
    _slider_simulationTimeline->setRange(0, (int)ceil(_solverSettings->timeEnd/_solverSettings->timeStepSize)); 
    _slider_simulationTimeline->setValue((int)ceil(_viewer->_simWorld->GetWorldTime()/_solverSettings->timeStepSize));
     _label_simulationTimeline = new QLabel(QString("%1").arg(_viewer->_simWorld->GetWorldTime())); 
     _label_simulationTimeline->setFixedWidth(200);
      _text_simulationTimeline = new QLabel("Simulation Time"); 
      _text_simulationTimeline->setFixedWidth(200);
    _button_resetSimulation = new QPushButton(this); 
    _button_resetSimulation->setText("Reset Simulation Time");
    _button_generateSlice_x = new QPushButton(this); 
    _button_generateSlice_x->setText("Generate x-slice");
    _button_generateSlice_y = new QPushButton(this); 
    _button_generateSlice_y->setText("Generate y-slice");
    _button_generateSlice_z = new QPushButton(this); 
    _button_generateSlice_z->setText("Generate z-slice");
    _button_clearPressures = new QPushButton(this); 
    _button_clearPressures->setText("Clear Pressure Field");
    _layout->addWidget(_viewer.get()      , 0, 0, 1, 4);
    _layout->addWidget(  _text_impulseScaling, 1, 0);
    _layout->addWidget(_slider_impulseScaling, 1, 1);
    _layout->addWidget( _label_impulseScaling, 1, 2);
    _layout->addWidget(  _text_simulationTimeline, 2, 0);
    _layout->addWidget(_slider_simulationTimeline, 2, 1);
    _layout->addWidget( _label_simulationTimeline, 2, 2);
    _layout->addWidget(_button_resetSimulation   , 2, 3);
    _layout->addLayout(_controlPanelLayout, 3, 0, 1, 4);
    _controlPanelLayout->addWidget(_button_generateSlice_x, 0, 0);
    _controlPanelLayout->addWidget(_button_generateSlice_y, 0, 1);
    _controlPanelLayout->addWidget(_button_generateSlice_z, 0, 2);
    _layout->addWidget(_button_clearPressures, 4, 0);
    setLayout(_layout);
    resize(800, 600);
    // signal-slot stuff
    connect(_slider_impulseScaling, SIGNAL(valueChanged(int)), this, SLOT(SliderValueChanged()));
    connect(_slider_simulationTimeline, SIGNAL(valueChanged(int)), this, SLOT(SliderValueChanged()));
    connect(_button_resetSimulation, SIGNAL(clicked()), this, SLOT(ResetSystemTime()));
    connect(_button_generateSlice_x, &QPushButton::clicked, this, [this]{
            GenerateSlice(0, _solverSettings->domainCenter[0]);});
    connect(_button_generateSlice_y, &QPushButton::clicked, this, [this]{
            GenerateSlice(1, _solverSettings->domainCenter[1]);});
    connect(_button_generateSlice_z, &QPushButton::clicked, this, [this]{
            GenerateSlice(2, _solverSettings->domainCenter[2]);});
    connect(_button_clearPressures, &QPushButton::clicked, this, [this]{
            ResetSystemTime(false);}); 
            
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
    if (_viewer->animationIsStarted())
        _viewer->stopAnimation(); 
    const REAL impulseScaling = (REAL)_slider_impulseScaling->value()*IMP_SLIDER_SCALE; 
    const REAL newTime = (REAL)_slider_simulationTimeline->value()*_solverSettings->timeStepSize; 
    _label_impulseScaling->setText(QString("%1").arg(impulseScaling));
    _label_simulationTimeline->setText(QString("%1").arg(newTime));
    // update simworld
    _viewer->_drawImpulseScaling = pow(10., impulseScaling);
    _viewer->_simWorld->SetWorldTime(newTime); 
    _viewer->_simWorld->UpdateObjectState(newTime); 
    _viewer->updateGL();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Widget::
ResetSystemTime(const bool &fromSlider)
{
    const REAL newTime = fromSlider ? 
        (REAL)_slider_simulationTimeline->value()*_solverSettings->timeStepSize
      : _viewer->_simWorld->GetWorldTime(); 
    _viewer->_simWorld->ResetStartTime(newTime);
    _viewer->SetAllSliceDataReady(false);
    _viewer->_currentFrame = 0;
    _viewer->updateGL();
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Widget::
GenerateSlice(const int &dim, const REAL &offset)
{
    _viewer->AddSlice(dim, offset);
    _viewer->updateGL();
}
