/********************************************************************************
** Form generated from reading UI file 'RigidSim.ui'
**
** Created by: Qt User Interface Compiler version 5.3.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RIGIDSIM_H
#define UI_RIGIDSIM_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QWidget>
#include "tools/rigidsim/RigidCanvas.h"

QT_BEGIN_NAMESPACE

class Ui_RigidSim
{
public:
    QAction *actionLoadBody;
    QAction *actionExportMatlab;
    QAction *actionExportAbaqus;
    QAction *actionModeShapes;
    QAction *actionLoadModes;
    QAction *actionExportMatrices;
    QAction *actionExportSurfaceMesh;
    QAction *actionExportImpulses;
    QAction *actionManual_Adjusting;
    QAction *action_Load;
    QAction *actionExit;
    QAction *actionStart;
    QAction *actionRecord_Positions;
    QAction *actionShakingPiggy;
    QAction *actionStep;
    QAction *actionDrop_Multiple_Objs;
    QAction *actionDropObjects;
    QAction *actionDropObjsWithFixed;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    RigidCanvas *canvas;
    QMenuBar *menubar;
    QMenu *menu_Demo;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *RigidSim)
    {
        if (RigidSim->objectName().isEmpty())
            RigidSim->setObjectName(QStringLiteral("RigidSim"));
        RigidSim->resize(889, 702);
        QIcon icon;
        icon.addFile(QStringLiteral(":/images/fracture.png"), QSize(), QIcon::Normal, QIcon::Off);
        RigidSim->setWindowIcon(icon);
        actionLoadBody = new QAction(RigidSim);
        actionLoadBody->setObjectName(QStringLiteral("actionLoadBody"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/images/load_frac.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoadBody->setIcon(icon1);
        actionExportMatlab = new QAction(RigidSim);
        actionExportMatlab->setObjectName(QStringLiteral("actionExportMatlab"));
        actionExportAbaqus = new QAction(RigidSim);
        actionExportAbaqus->setObjectName(QStringLiteral("actionExportAbaqus"));
        actionModeShapes = new QAction(RigidSim);
        actionModeShapes->setObjectName(QStringLiteral("actionModeShapes"));
        actionLoadModes = new QAction(RigidSim);
        actionLoadModes->setObjectName(QStringLiteral("actionLoadModes"));
        actionExportMatrices = new QAction(RigidSim);
        actionExportMatrices->setObjectName(QStringLiteral("actionExportMatrices"));
        actionExportSurfaceMesh = new QAction(RigidSim);
        actionExportSurfaceMesh->setObjectName(QStringLiteral("actionExportSurfaceMesh"));
        actionExportImpulses = new QAction(RigidSim);
        actionExportImpulses->setObjectName(QStringLiteral("actionExportImpulses"));
        actionManual_Adjusting = new QAction(RigidSim);
        actionManual_Adjusting->setObjectName(QStringLiteral("actionManual_Adjusting"));
        actionManual_Adjusting->setCheckable(true);
        action_Load = new QAction(RigidSim);
        action_Load->setObjectName(QStringLiteral("action_Load"));
        actionExit = new QAction(RigidSim);
        actionExit->setObjectName(QStringLiteral("actionExit"));
        actionStart = new QAction(RigidSim);
        actionStart->setObjectName(QStringLiteral("actionStart"));
        actionRecord_Positions = new QAction(RigidSim);
        actionRecord_Positions->setObjectName(QStringLiteral("actionRecord_Positions"));
        actionRecord_Positions->setCheckable(true);
        actionShakingPiggy = new QAction(RigidSim);
        actionShakingPiggy->setObjectName(QStringLiteral("actionShakingPiggy"));
        actionStep = new QAction(RigidSim);
        actionStep->setObjectName(QStringLiteral("actionStep"));
        actionDrop_Multiple_Objs = new QAction(RigidSim);
        actionDrop_Multiple_Objs->setObjectName(QStringLiteral("actionDrop_Multiple_Objs"));
        actionDropObjects = new QAction(RigidSim);
        actionDropObjects->setObjectName(QStringLiteral("actionDropObjects"));
        actionDropObjsWithFixed = new QAction(RigidSim);
        actionDropObjsWithFixed->setObjectName(QStringLiteral("actionDropObjsWithFixed"));
        centralwidget = new QWidget(RigidSim);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        horizontalLayout = new QHBoxLayout(centralwidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        canvas = new RigidCanvas(centralwidget);
        canvas->setObjectName(QStringLiteral("canvas"));

        horizontalLayout->addWidget(canvas);

        RigidSim->setCentralWidget(centralwidget);
        menubar = new QMenuBar(RigidSim);
        menubar->setObjectName(QStringLiteral("menubar"));
        menubar->setGeometry(QRect(0, 0, 889, 26));
        menu_Demo = new QMenu(menubar);
        menu_Demo->setObjectName(QStringLiteral("menu_Demo"));
        RigidSim->setMenuBar(menubar);
        statusbar = new QStatusBar(RigidSim);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        RigidSim->setStatusBar(statusbar);

        menubar->addAction(menu_Demo->menuAction());
        menu_Demo->addAction(actionDropObjects);
        menu_Demo->addAction(actionDropObjsWithFixed);
        menu_Demo->addAction(actionShakingPiggy);
        menu_Demo->addSeparator();
        menu_Demo->addAction(actionStart);
        menu_Demo->addAction(actionStep);
        menu_Demo->addSeparator();
        menu_Demo->addAction(actionExit);

        retranslateUi(RigidSim);

        QMetaObject::connectSlotsByName(RigidSim);
    } // setupUi

    void retranslateUi(QMainWindow *RigidSim)
    {
        RigidSim->setWindowTitle(QApplication::translate("RigidSim", "Rigid Simulator", 0));
        actionLoadBody->setText(QApplication::translate("RigidSim", "Load Rigid Body", 0));
        actionExportMatlab->setText(QApplication::translate("RigidSim", "&Export Matlab Data", 0));
        actionExportAbaqus->setText(QApplication::translate("RigidSim", "&Export Abaqus Mesh", 0));
        actionModeShapes->setText(QApplication::translate("RigidSim", "Mode Shapes", 0));
        actionLoadModes->setText(QApplication::translate("RigidSim", "Load Modes", 0));
        actionExportMatrices->setText(QApplication::translate("RigidSim", "Export Matrices", 0));
        actionExportSurfaceMesh->setText(QApplication::translate("RigidSim", "Export Surface Mesh", 0));
        actionExportImpulses->setText(QApplication::translate("RigidSim", "Export Impulses", 0));
        actionManual_Adjusting->setText(QApplication::translate("RigidSim", "Manual Adjusting", 0));
        action_Load->setText(QApplication::translate("RigidSim", "&Load", 0));
        actionExit->setText(QApplication::translate("RigidSim", "&Exit", 0));
        actionStart->setText(QApplication::translate("RigidSim", "&Start", 0));
        actionStart->setShortcut(QApplication::translate("RigidSim", "F8", 0));
        actionRecord_Positions->setText(QApplication::translate("RigidSim", "Record Positions", 0));
        actionShakingPiggy->setText(QApplication::translate("RigidSim", "Shaking piggy bank", 0));
        actionStep->setText(QApplication::translate("RigidSim", "Step", 0));
        actionStep->setShortcut(QApplication::translate("RigidSim", "F5", 0));
        actionDrop_Multiple_Objs->setText(QApplication::translate("RigidSim", "Drop Multiple Objs", 0));
        actionDropObjects->setText(QApplication::translate("RigidSim", "Drop objects", 0));
        actionDropObjsWithFixed->setText(QApplication::translate("RigidSim", "Drop objs with fixed ", 0));
        menu_Demo->setTitle(QApplication::translate("RigidSim", "&Demo", 0));
    } // retranslateUi

};

namespace Ui {
    class RigidSim: public Ui_RigidSim {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RIGIDSIM_H
