#include <fstream> 
#include <iostream> 
#include <Eigen/Dense> 
#include <cstdlib> 
#include <ctime> 

#include "resample.h" 
#include "Texture.h" 

#include <QApplication> 
#include <QMainWindow> 
#include "qcustomplot.h" 

#include "convert.h" 


using namespace std; 
using namespace SIGNAL_PROCESSING;

int main(int argc, char *argv[])
{




    SoundTexture texture( "../SADDES_smoothBL4_v20mph_y90_finalSound_nonextended.wav" ); 

    Eigen::VectorXd L = texture.L(); 
    Eigen::VectorXd R = texture.R(); 

    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(L.size(), 0, 0.0001*L.size());


    Eigen::MatrixXd tmp = texture.BuildGraph();

    cout << tmp << endl;

    QApplication a(argc, argv); 
    QMainWindow window; 

    QCustomPlot customPlot; 
    window.setCentralWidget( &customPlot ); 

    QCPPlotTitle *title = new QCPPlotTitle( &customPlot ); 
    title->setText( "sound signal" );
    title->setFont(QFont("sans", 12, QFont::Bold));

    customPlot.plotLayout()->addElement(0, 0, NULL);
    customPlot.addGraph(); 
    customPlot.graph(0)->setData( CONVERT::Eigenvec2QVec(t), CONVERT::Eigenvec2QVec(L) ); 

    customPlot.xAxis->setLabel("x");
    customPlot.yAxis->setLabel("y");
    customPlot.rescaleAxes();

    window.setGeometry(5000, 500, 500, 400);
    window.show();
    return a.exec();

}
