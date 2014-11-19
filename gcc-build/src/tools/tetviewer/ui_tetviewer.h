/********************************************************************************
** Form generated from reading UI file 'tetviewer.ui'
**
** Created by: Qt User Interface Compiler version 5.3.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TETVIEWER_H
#define UI_TETVIEWER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDial>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "tools/tetviewer/TetViewerCanvas.h"

QT_BEGIN_NAMESPACE

class Ui_TetViewerFrame
{
public:
    QAction *actionOpen;
    QAction *actionModeShapes;
    QAction *actionLoadModes;
    QAction *actionWireframe;
    QAction *actionBinaryTetFormat;
    QAction *actionTextTetFormat;
    QAction *actionMeshInfo;
    QAction *actionCheckUselessVertex;
    QAction *actionAbaqusTetFormat;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    QWidget *buttonContainer;
    QVBoxLayout *verticalLayout;
    QLabel *modeIndexLabel;
    QSpinBox *modeIndex;
    QSpacerItem *horizontalSpacer;
    QFrame *line;
    QLabel *modeCoordinateLabel;
    QDial *modalCoordinate;
    QLabel *modeScaleLabel;
    QDoubleSpinBox *modeScale;
    QSpacerItem *horizontalSpacer_2;
    QFrame *line_2;
    QLabel *objectDensityLabel;
    QDoubleSpinBox *objectDensity;
    TetViewerCanvas *canvas;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuExport;
    QMenu *menuModal;
    QMenu *menuView;
    QMenu *menuMesh;
    QStatusBar *statusbar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *TetViewerFrame)
    {
        if (TetViewerFrame->objectName().isEmpty())
            TetViewerFrame->setObjectName(QStringLiteral("TetViewerFrame"));
        TetViewerFrame->resize(1280, 702);
        QIcon icon;
        icon.addFile(QStringLiteral(":/images/obb.png"), QSize(), QIcon::Normal, QIcon::Off);
        TetViewerFrame->setWindowIcon(icon);
        actionOpen = new QAction(TetViewerFrame);
        actionOpen->setObjectName(QStringLiteral("actionOpen"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/images/load.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpen->setIcon(icon1);
        actionModeShapes = new QAction(TetViewerFrame);
        actionModeShapes->setObjectName(QStringLiteral("actionModeShapes"));
        actionLoadModes = new QAction(TetViewerFrame);
        actionLoadModes->setObjectName(QStringLiteral("actionLoadModes"));
        QIcon icon2;
        icon2.addFile(QStringLiteral(":/images/load_frac.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoadModes->setIcon(icon2);
        actionWireframe = new QAction(TetViewerFrame);
        actionWireframe->setObjectName(QStringLiteral("actionWireframe"));
        actionWireframe->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QStringLiteral(":/images/wireframe.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionWireframe->setIcon(icon3);
        actionBinaryTetFormat = new QAction(TetViewerFrame);
        actionBinaryTetFormat->setObjectName(QStringLiteral("actionBinaryTetFormat"));
        QIcon icon4;
        icon4.addFile(QStringLiteral(":/images/binary.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionBinaryTetFormat->setIcon(icon4);
        actionTextTetFormat = new QAction(TetViewerFrame);
        actionTextTetFormat->setObjectName(QStringLiteral("actionTextTetFormat"));
        actionMeshInfo = new QAction(TetViewerFrame);
        actionMeshInfo->setObjectName(QStringLiteral("actionMeshInfo"));
        actionMeshInfo->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QStringLiteral(":/images/info.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMeshInfo->setIcon(icon5);
        actionCheckUselessVertex = new QAction(TetViewerFrame);
        actionCheckUselessVertex->setObjectName(QStringLiteral("actionCheckUselessVertex"));
        actionAbaqusTetFormat = new QAction(TetViewerFrame);
        actionAbaqusTetFormat->setObjectName(QStringLiteral("actionAbaqusTetFormat"));
        centralwidget = new QWidget(TetViewerFrame);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        horizontalLayout = new QHBoxLayout(centralwidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMinAndMaxSize);
        buttonContainer = new QWidget(centralwidget);
        buttonContainer->setObjectName(QStringLiteral("buttonContainer"));
        buttonContainer->setEnabled(true);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(buttonContainer->sizePolicy().hasHeightForWidth());
        buttonContainer->setSizePolicy(sizePolicy);
        buttonContainer->setMinimumSize(QSize(300, 0));
        buttonContainer->setMaximumSize(QSize(300, 16777215));
        buttonContainer->setBaseSize(QSize(300, 597));
        verticalLayout = new QVBoxLayout(buttonContainer);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetMinAndMaxSize);
        modeIndexLabel = new QLabel(buttonContainer);
        modeIndexLabel->setObjectName(QStringLiteral("modeIndexLabel"));
        modeIndexLabel->setEnabled(true);

        verticalLayout->addWidget(modeIndexLabel);

        modeIndex = new QSpinBox(buttonContainer);
        modeIndex->setObjectName(QStringLiteral("modeIndex"));
        modeIndex->setEnabled(false);

        verticalLayout->addWidget(modeIndex);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        verticalLayout->addItem(horizontalSpacer);

        line = new QFrame(buttonContainer);
        line->setObjectName(QStringLiteral("line"));
        line->setLineWidth(10);
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        modeCoordinateLabel = new QLabel(buttonContainer);
        modeCoordinateLabel->setObjectName(QStringLiteral("modeCoordinateLabel"));
        modeCoordinateLabel->setTextFormat(Qt::AutoText);

        verticalLayout->addWidget(modeCoordinateLabel);

        modalCoordinate = new QDial(buttonContainer);
        modalCoordinate->setObjectName(QStringLiteral("modalCoordinate"));
        modalCoordinate->setEnabled(false);
        modalCoordinate->setMinimum(-100);
        modalCoordinate->setMaximum(100);

        verticalLayout->addWidget(modalCoordinate);

        modeScaleLabel = new QLabel(buttonContainer);
        modeScaleLabel->setObjectName(QStringLiteral("modeScaleLabel"));

        verticalLayout->addWidget(modeScaleLabel);

        modeScale = new QDoubleSpinBox(buttonContainer);
        modeScale->setObjectName(QStringLiteral("modeScale"));
        modeScale->setEnabled(false);
        modeScale->setDecimals(6);
        modeScale->setMaximum(1000);
        modeScale->setSingleStep(0.1);
        modeScale->setValue(1);

        verticalLayout->addWidget(modeScale);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        verticalLayout->addItem(horizontalSpacer_2);

        line_2 = new QFrame(buttonContainer);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setLineWidth(10);
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line_2);

        objectDensityLabel = new QLabel(buttonContainer);
        objectDensityLabel->setObjectName(QStringLiteral("objectDensityLabel"));

        verticalLayout->addWidget(objectDensityLabel);

        objectDensity = new QDoubleSpinBox(buttonContainer);
        objectDensity->setObjectName(QStringLiteral("objectDensity"));
        objectDensity->setMaximum(20000);
        objectDensity->setValue(1000);

        verticalLayout->addWidget(objectDensity);


        horizontalLayout->addWidget(buttonContainer);

        canvas = new TetViewerCanvas(centralwidget);
        canvas->setObjectName(QStringLiteral("canvas"));
        canvas->setMinimumSize(QSize(956, 0));

        horizontalLayout->addWidget(canvas);

        TetViewerFrame->setCentralWidget(centralwidget);
        menubar = new QMenuBar(TetViewerFrame);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 1280, 25));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QStringLiteral("menuFile"));
        menuExport = new QMenu(menuFile);
        menuExport->setObjectName(QStringLiteral("menuExport"));
        menuModal = new QMenu(menubar);
        menuModal->setObjectName(QStringLiteral("menuModal"));
        menuView = new QMenu(menubar);
        menuView->setObjectName(QStringLiteral("menuView"));
        menuMesh = new QMenu(menubar);
        menuMesh->setObjectName(QStringLiteral("menuMesh"));
        TetViewerFrame->setMenuBar(menubar);
        statusbar = new QStatusBar(TetViewerFrame);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        TetViewerFrame->setStatusBar(statusbar);
        toolBar = new QToolBar(TetViewerFrame);
        toolBar->setObjectName(QStringLiteral("toolBar"));
        TetViewerFrame->addToolBar(Qt::TopToolBarArea, toolBar);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuView->menuAction());
        menubar->addAction(menuMesh->menuAction());
        menubar->addAction(menuModal->menuAction());
        menuFile->addAction(actionOpen);
        menuFile->addAction(menuExport->menuAction());
        menuExport->addAction(actionBinaryTetFormat);
        menuExport->addAction(actionTextTetFormat);
        menuExport->addAction(actionAbaqusTetFormat);
        menuModal->addSeparator();
        menuModal->addAction(actionLoadModes);
        menuModal->addAction(actionModeShapes);
        menuView->addAction(actionWireframe);
        menuView->addSeparator();
        menuView->addAction(actionMeshInfo);
        menuMesh->addAction(actionCheckUselessVertex);
        toolBar->addAction(actionOpen);
        toolBar->addAction(actionLoadModes);
        toolBar->addAction(actionBinaryTetFormat);
        toolBar->addSeparator();
        toolBar->addAction(actionWireframe);
        toolBar->addAction(actionMeshInfo);

        retranslateUi(TetViewerFrame);

        QMetaObject::connectSlotsByName(TetViewerFrame);
    } // setupUi

    void retranslateUi(QMainWindow *TetViewerFrame)
    {
        TetViewerFrame->setWindowTitle(QApplication::translate("TetViewerFrame", "Tetrahedron Viewer", 0));
        actionOpen->setText(QApplication::translate("TetViewerFrame", "&Open", 0));
        actionModeShapes->setText(QApplication::translate("TetViewerFrame", "Mode Shapes", 0));
        actionLoadModes->setText(QApplication::translate("TetViewerFrame", "Load Modes", 0));
        actionWireframe->setText(QApplication::translate("TetViewerFrame", "&Wireframe", 0));
        actionBinaryTetFormat->setText(QApplication::translate("TetViewerFrame", "Binary Tet Format (.tet)", 0));
        actionTextTetFormat->setText(QApplication::translate("TetViewerFrame", "Text Tet Format (.node/.ele)", 0));
        actionMeshInfo->setText(QApplication::translate("TetViewerFrame", "&Mesh Info", 0));
        actionCheckUselessVertex->setText(QApplication::translate("TetViewerFrame", "Check Useless Vertex", 0));
        actionAbaqusTetFormat->setText(QApplication::translate("TetViewerFrame", "Abaqus Tet Format (.aba)", 0));
        modeIndexLabel->setText(QApplication::translate("TetViewerFrame", "Mode to display", 0));
        modeCoordinateLabel->setText(QApplication::translate("TetViewerFrame", "Modal coordinate", 0));
        modeScaleLabel->setText(QApplication::translate("TetViewerFrame", "Mode scale", 0));
        objectDensityLabel->setText(QApplication::translate("TetViewerFrame", "Object Density (kg/m^3)", 0));
        menuFile->setTitle(QApplication::translate("TetViewerFrame", "&File", 0));
        menuExport->setTitle(QApplication::translate("TetViewerFrame", "Export", 0));
        menuModal->setTitle(QApplication::translate("TetViewerFrame", "&Modal", 0));
        menuView->setTitle(QApplication::translate("TetViewerFrame", "&View", 0));
        menuMesh->setTitle(QApplication::translate("TetViewerFrame", "M&esh", 0));
        toolBar->setWindowTitle(QApplication::translate("TetViewerFrame", "toolBar", 0));
    } // retranslateUi

};

namespace Ui {
    class TetViewerFrame: public Ui_TetViewerFrame {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TETVIEWER_H
