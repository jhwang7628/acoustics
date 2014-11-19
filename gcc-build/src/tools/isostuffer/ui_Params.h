/********************************************************************************
** Form generated from reading UI file 'Params.ui'
**
** Created by: Qt User Interface Compiler version 5.3.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PARAMS_H
#define UI_PARAMS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_ParamsDialog
{
public:
    QVBoxLayout *verticalLayout_2;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_3;
    QSpinBox *resSpin;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_4;
    QSpinBox *levelSpin;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_5;
    QSpinBox *marginSpin;
    QHBoxLayout *horizontalLayout_2;
    QLabel *alphaLongLabel;
    QLineEdit *alphaLongBox;
    QHBoxLayout *horizontalLayout;
    QLabel *alphaShortLabel;
    QLineEdit *alphaShortBox;

    void setupUi(QDialog *ParamsDialog)
    {
        if (ParamsDialog->objectName().isEmpty())
            ParamsDialog->setObjectName(QStringLiteral("ParamsDialog"));
        ParamsDialog->resize(287, 220);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ParamsDialog->sizePolicy().hasHeightForWidth());
        ParamsDialog->setSizePolicy(sizePolicy);
        ParamsDialog->setMinimumSize(QSize(287, 190));
        ParamsDialog->setMaximumSize(QSize(287, 220));
        verticalLayout_2 = new QVBoxLayout(ParamsDialog);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        groupBox = new QGroupBox(ParamsDialog);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_4->addWidget(label_3);

        resSpin = new QSpinBox(groupBox);
        resSpin->setObjectName(QStringLiteral("resSpin"));
        resSpin->setMinimum(2);
        resSpin->setMaximum(8);
        resSpin->setValue(5);

        horizontalLayout_4->addWidget(resSpin);


        verticalLayout->addLayout(horizontalLayout_4);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_3->addWidget(label_4);

        levelSpin = new QSpinBox(groupBox);
        levelSpin->setObjectName(QStringLiteral("levelSpin"));
        levelSpin->setLayoutDirection(Qt::LeftToRight);
        levelSpin->setMinimum(1);
        levelSpin->setMaximum(5);
        levelSpin->setValue(3);

        horizontalLayout_3->addWidget(levelSpin);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_5 = new QLabel(groupBox);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setLayoutDirection(Qt::LeftToRight);
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        horizontalLayout_5->addWidget(label_5);

        marginSpin = new QSpinBox(groupBox);
        marginSpin->setObjectName(QStringLiteral("marginSpin"));
        marginSpin->setLayoutDirection(Qt::LeftToRight);
        marginSpin->setMinimum(1);
        marginSpin->setMaximum(100);
        marginSpin->setValue(7);

        horizontalLayout_5->addWidget(marginSpin);


        verticalLayout->addLayout(horizontalLayout_5);


        verticalLayout_2->addWidget(groupBox);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        alphaLongLabel = new QLabel(ParamsDialog);
        alphaLongLabel->setObjectName(QStringLiteral("alphaLongLabel"));

        horizontalLayout_2->addWidget(alphaLongLabel);

        alphaLongBox = new QLineEdit(ParamsDialog);
        alphaLongBox->setObjectName(QStringLiteral("alphaLongBox"));
        alphaLongBox->setLayoutDirection(Qt::RightToLeft);

        horizontalLayout_2->addWidget(alphaLongBox);


        verticalLayout_2->addLayout(horizontalLayout_2);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        alphaShortLabel = new QLabel(ParamsDialog);
        alphaShortLabel->setObjectName(QStringLiteral("alphaShortLabel"));

        horizontalLayout->addWidget(alphaShortLabel);

        alphaShortBox = new QLineEdit(ParamsDialog);
        alphaShortBox->setObjectName(QStringLiteral("alphaShortBox"));
        alphaShortBox->setLayoutDirection(Qt::RightToLeft);

        horizontalLayout->addWidget(alphaShortBox);


        verticalLayout_2->addLayout(horizontalLayout);


        retranslateUi(ParamsDialog);

        QMetaObject::connectSlotsByName(ParamsDialog);
    } // setupUi

    void retranslateUi(QDialog *ParamsDialog)
    {
        ParamsDialog->setWindowTitle(QApplication::translate("ParamsDialog", "Parameters", 0));
        groupBox->setTitle(QApplication::translate("ParamsDialog", "Oct Tree", 0));
        label_3->setText(QApplication::translate("ParamsDialog", "Resolution:  ", 0));
        label_4->setText(QApplication::translate("ParamsDialog", "# of Level:  ", 0));
        label_5->setText(QApplication::translate("ParamsDialog", "margin:  ", 0));
        alphaLongLabel->setText(QApplication::translate("ParamsDialog", "  Alpha (long)  : ", 0));
        alphaLongBox->setText(QApplication::translate("ParamsDialog", "0.25", 0));
        alphaShortLabel->setText(QApplication::translate("ParamsDialog", "  Alpha (short): ", 0));
        alphaShortBox->setText(QApplication::translate("ParamsDialog", "0.42978", 0));
    } // retranslateUi

};

namespace Ui {
    class ParamsDialog: public Ui_ParamsDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PARAMS_H
