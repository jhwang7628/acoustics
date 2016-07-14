#include <math.h>
#include <stdlib.h> // RAND_MAX
#include <ui/ModalViewer.h>

using namespace qglviewer;
using namespace std;

//##############################################################################
//##############################################################################
void ModalViewer::init()
{
    restoreStateFromFile();
    glDisable(GL_LIGHTING);
    glPointSize(3.0);
    setGridIsDrawn();
    help();
    startAnimation();
}

//##############################################################################
//##############################################################################
void ModalViewer::draw()
{
    glBegin(GL_POINTS);
    glEnd();
}

//##############################################################################
//##############################################################################
void ModalViewer::animate()
{
}

//##############################################################################
//##############################################################################
QString ModalViewer::helpString() const
{
    QString text("<h2>Modal ModalViewer</h2>");
    text += "Used for debugging FDTD_RigidSoundObject and associated classes";
    return text;
}

