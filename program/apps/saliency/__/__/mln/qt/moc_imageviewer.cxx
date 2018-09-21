/****************************************************************************
** Meta object code from reading C++ file 'imageviewer.hpp'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../../../mln/qt/imageviewer.hpp"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'imageviewer.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_mln__qt__ImageViewer[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      25,   22,   21,   21, 0x05,
      53,   22,   21,   21, 0x05,

 // slots: signature, parameters, type, tag, flags
      73,   21,   21,   21, 0x0a,
      82,   21,   21,   21, 0x0a,
      96,   91,   21,   21, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_mln__qt__ImageViewer[] = {
    "mln::qt::ImageViewer\0\0pt\0"
    "pointSelected(mln::point2d)\0"
    "pointHover(point2d)\0update()\0notify()\0"
    "rect\0onzoom(QRectF)\0"
};

void mln::qt::ImageViewer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ImageViewer *_t = static_cast<ImageViewer *>(_o);
        switch (_id) {
        case 0: _t->pointSelected((*reinterpret_cast< const mln::point2d(*)>(_a[1]))); break;
        case 1: _t->pointHover((*reinterpret_cast< const point2d(*)>(_a[1]))); break;
        case 2: _t->update(); break;
        case 3: _t->notify(); break;
        case 4: _t->onzoom((*reinterpret_cast< const QRectF(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData mln::qt::ImageViewer::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject mln::qt::ImageViewer::staticMetaObject = {
    { &QGraphicsView::staticMetaObject, qt_meta_stringdata_mln__qt__ImageViewer,
      qt_meta_data_mln__qt__ImageViewer, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &mln::qt::ImageViewer::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *mln::qt::ImageViewer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *mln::qt::ImageViewer::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_mln__qt__ImageViewer))
        return static_cast<void*>(const_cast< ImageViewer*>(this));
    return QGraphicsView::qt_metacast(_clname);
}

int mln::qt::ImageViewer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGraphicsView::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    }
    return _id;
}

// SIGNAL 0
void mln::qt::ImageViewer::pointSelected(const mln::point2d & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void mln::qt::ImageViewer::pointHover(const point2d & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
