#include "brush.hpp"
#include "constants.hpp"
#include <QtGui>
#include <iostream>
#include <mln/io/imsave.hpp>

inline
void brush(QPixmap* pixmap, const QPointF& position, QColor color, int r)
{
  QPainter painter(pixmap);
  painter.setBrush(color);
  painter.setPen(Qt::NoPen);

  QRectF rectangle(position.x() - r, position.y() - r, 2*r, 2*r);
  painter.drawEllipse(rectangle);
}



MyBrush::MyBrush(mln::qt::ImageViewer* viewer,
                 callback_fun_t callback)
    : m_viewer (viewer),
      m_callback(callback),
      m_scene(viewer->getScene()),
      m_pixmap(viewer->getPixmap()),
      m_ima(m_pixmap->pixmap()),
      m_active(false),
      m_radius(5),
      m_reject_value (0.5),
      m_policy (0)
{
  m_ori = mln::clone(viewer->getView());
}

bool
MyBrush::eventFilter(QObject* obj, QEvent* ev)
{
  (void) obj;
  if (m_active and ev->type() == QEvent::GraphicsSceneMouseMove)
    {
      QGraphicsSceneMouseEvent* event = static_cast<QGraphicsSceneMouseEvent*>(ev);
      if (event->buttons() & Qt::LeftButton)
        {
          brush(&m_ima, m_pixmap->mapFromScene(event->scenePos()), m_color, m_radius);
          m_pixmap->setPixmap(m_ima);
        }
    }
  return false;
}

void
MyBrush::setColor1()
{
  m_active = true;
  m_color = Qt::red;
}

void
MyBrush::setColor2()
{
  m_active = true;
  m_color = Qt::blue;
}


void
MyBrush::setRadius(int k)
{
  m_radius = k;
}


void
MyBrush::run()
{
  std::cout << "Running." << std::endl;
  std::cout << "Reject value: " << m_reject_value << std::endl;
  std::cout << "Policy: " << std::hex << m_policy << std::endl;
  m_viewer->notify();


  mln::image2d<mln::rgb8>& view = m_viewer->getView();
  view = m_callback(view, m_reject_value, m_policy);
  m_viewer->update();
}

void
MyBrush::revert()
{
  m_pixmap->setPixmap(m_ima);
}


void
MyBrush::reload()
{
  m_viewer->getView() = mln::clone(m_ori);
  m_viewer->update();
  m_ima = m_pixmap->pixmap(); // FLUSH THE MARKERS
}

void
MyBrush::set_reject_value(const QString& v)
{
  m_reject_value = v.toFloat();
}

void
MyBrush::set_ambiguity_policy(int index)
{
  QComboBox* combo = qobject_cast< QComboBox * >(sender());
  if (combo == NULL)
    return;

  QVariant data = combo->itemData(index);
  m_policy = m_policy & ~AMBIGUITY_MASK;
  m_policy |= data.toUInt();
}
