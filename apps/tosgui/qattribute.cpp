#include "qattribute.hpp"
#include <QEvent>
#include <QMouseEvent>
#include <qwt_picker_machine.h>

namespace mln
{


  QAttributeBase::QAttributeBase(const QwtText& name)
    : QwtPlot(name)
  {
    m_curve = new QwtPlotCurve("Energy");
    m_curve->setSamples(m_data);
    m_curve->setStyle(QwtPlotCurve::Lines);
    m_curve->attach(this);

    std::cout << this->canvas() << std::endl;

    picker = new QwtPlotPicker(this->canvas());
    picker->setStateMachine(new  QwtPickerTrackerMachine);
    picker->setTrackerMode(QwtPicker::AlwaysOn);
    picker->setRubberBand(QwtPicker::VLineRubberBand);
    picker->setEnabled(true);
    this->canvas()->installEventFilter(this);
  }

  void
  QAttributeBase::keyReleaseEvent(QKeyEvent* event)
  {
    if (event->key() == Qt::Key_F)
      this->showFilteringWindow();
  }

}
