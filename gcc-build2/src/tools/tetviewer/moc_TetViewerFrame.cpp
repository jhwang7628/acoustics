/****************************************************************************
** Meta object code from reading C++ file 'TetViewerFrame.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.3.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../src/tools/tetviewer/TetViewerFrame.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TetViewerFrame.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.3.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_TetViewerFrame_t {
    QByteArrayData data[9];
    char stringdata[127];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_TetViewerFrame_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_TetViewerFrame_t qt_meta_stringdata_TetViewerFrame = {
    {
QT_MOC_LITERAL(0, 0, 14),
QT_MOC_LITERAL(1, 15, 4),
QT_MOC_LITERAL(2, 20, 0),
QT_MOC_LITERAL(3, 21, 10),
QT_MOC_LITERAL(4, 32, 14),
QT_MOC_LITERAL(5, 47, 17),
QT_MOC_LITERAL(6, 65, 17),
QT_MOC_LITERAL(7, 83, 18),
QT_MOC_LITERAL(8, 102, 24)
    },
    "TetViewerFrame\0open\0\0load_modes\0"
    "export_bin_tet\0export_abaqus_tet\0"
    "check_useless_vtx\0update_active_mode\0"
    "update_mode_displacement"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_TetViewerFrame[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x0a /* Public */,
       3,    0,   50,    2, 0x0a /* Public */,
       4,    0,   51,    2, 0x0a /* Public */,
       5,    0,   52,    2, 0x0a /* Public */,
       6,    0,   53,    2, 0x0a /* Public */,
       7,    0,   54,    2, 0x0a /* Public */,
       8,    0,   55,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void TetViewerFrame::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        TetViewerFrame *_t = static_cast<TetViewerFrame *>(_o);
        switch (_id) {
        case 0: _t->open(); break;
        case 1: _t->load_modes(); break;
        case 2: _t->export_bin_tet(); break;
        case 3: _t->export_abaqus_tet(); break;
        case 4: _t->check_useless_vtx(); break;
        case 5: _t->update_active_mode(); break;
        case 6: _t->update_mode_displacement(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject TetViewerFrame::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_TetViewerFrame.data,
      qt_meta_data_TetViewerFrame,  qt_static_metacall, 0, 0}
};


const QMetaObject *TetViewerFrame::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *TetViewerFrame::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_TetViewerFrame.stringdata))
        return static_cast<void*>(const_cast< TetViewerFrame*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int TetViewerFrame::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
