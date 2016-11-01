#include <ui/FDTD_AcousticSimulator_Widget.h>

//##############################################################################
//##############################################################################
FDTD_AcousticSimulator_Widget::
FDTD_AcousticSimulator_Widget(std::shared_ptr<FDTD_AcousticSimulator_Viewer> &viewer)
    : _viewer(viewer)
{
    // initialize layout
    _layout                = new QGridLayout(); 
    _slider_impulseScaling = new QSlider(Qt::Horizontal, this); 
     _label_impulseScaling = new QLabel("0"); 
    _slider_impulseScaling->setRange(-3./IMP_SLIDER_SCALE, 2./IMP_SLIDER_SCALE); 
    _slider_impulseScaling->setValue(0);
    _layout->addWidget(_viewer.get()         , 0, 0, 1, 2);
    _layout->addWidget(_slider_impulseScaling, 1, 0);
    _layout->addWidget( _label_impulseScaling, 1, 1);
    setLayout(_layout);
    resize(800, 600);
    // signal-slot stuff
    connect(_slider_impulseScaling, SIGNAL(valueChanged(int)), this, SLOT(SliderValueChanged()));
}

//##############################################################################
//##############################################################################
void FDTD_AcousticSimulator_Widget::
SliderValueChanged()
{
    const REAL impulseScaling = (REAL)_slider_impulseScaling->value()*IMP_SLIDER_SCALE; 
    _label_impulseScaling->setText(QString("%1").arg(impulseScaling));
    _viewer->_drawImpulseScaling = pow(10., impulseScaling); 

    _viewer->updateGL();
}
