#include <mln/qt/imageviewer.hpp>
#include <iostream>
#include <QtGui>

namespace mln
{

  namespace qt
  {

    class ImageViewerEventHandler : public QObject
    {
      Q_OBJECT;
    public:
      ImageViewerEventHandler(ImageViewer* viewer)
      : m_viewer(viewer),
        m_view_center(0,0)
      {
      }

      void zoomBy(float factor);
      void startZoomRect(QPointF wpos, QPointF scene_pos);
      void endZoomRect(QPointF wpos, QPointF scene_pos);

    protected:
      virtual bool eventFilter(QObject* obj, QEvent* event);

    private:
      ImageViewer*      m_viewer;

      // Structure that stores zooming information
      struct {
        bool            active = false;
        QPointF         selection_start;  // Scene position
        QPointF         selection_wstart; // Widget position
      } zooming;
      QPointF           m_view_center;          // center of the viewport
      bool              m_auto_resize_window;   // resize the window when zooming
    };

    void
    ImageViewerEventHandler::zoomBy(float factor)
    {
      float scale = m_viewer->m_pixmap.scale() * factor;
      if (scale < 0.01) scale = 0.01;
      else if (scale > 100) scale = 100;
      m_view_center *= (scale /  m_viewer->m_pixmap.scale());

      // if (m_auto_resize_window)
      //   m_viewer->viewport()->
      // else
      //   m_viewer->viewport()->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

      m_viewer->m_pixmap.setScale(scale);
      m_viewer->m_scene->setSceneRect(m_viewer->m_pixmap.sceneBoundingRect());
      m_viewer->centerOn(m_view_center);
    }

    void
    ImageViewerEventHandler::startZoomRect(QPointF wpos, QPointF scene_pos)
    {
      zooming.active = true;
      zooming.selection_start = scene_pos;
      zooming.selection_wstart = wpos;
    }

    void
    ImageViewerEventHandler::endZoomRect(QPointF, QPointF scene_pos)
    {
      zooming.active = false;

      QPointF p1 = zooming.selection_start;
      QPointF p2 = scene_pos;
      int w = m_viewer->width();
      int h = m_viewer->height();
      qreal x1, x2, y1, y2;
      std::tie(x1,x2) = std::minmax(p1.x(), p2.x());
      std::tie(y1,y2) = std::minmax(p1.y(), p2.y());
      QPoint d = m_viewer->mapFromScene(x2-x1, y2-y1);
      float wscale = w / d.x();
      float yscale = h / d.y();
      float ratio = float(d.y()) / d.x();
      m_view_center.setX(0.5 * (x1+x2));
      m_view_center.setY(0.5 * (y1+y2));
      this->zoomBy(yscale);
      m_viewer->resize(w, w*ratio);
      m_viewer->centerOn(m_view_center);
      //m_viewer->centerOn(0.5 * yscale * (x1+x2), 0.5 * yscale * (y1+y2));
    }


    bool
    ImageViewerEventHandler::eventFilter(QObject* obj, QEvent* ev)
    {
      if (ev->type() == QEvent::KeyPress)
        {
          QKeyEvent* event = static_cast<QKeyEvent*>(ev);
          switch (event->key())
            {
              case '<': m_auto_resize_window = true; zoomBy(0.5); return true;
              case '>': m_auto_resize_window = true; zoomBy(2.0); return true;
            }
        }
       else if (ev->type() == QEvent::GraphicsSceneMousePress)
         {
           QGraphicsSceneMouseEvent* event = static_cast<QGraphicsSceneMouseEvent*>(ev);
           if (event->button() == Qt::RightButton)
             this->startZoomRect(event->pos(), event->scenePos());
         }
       else if (ev->type() == QEvent::GraphicsSceneMouseRelease)
         {
           QGraphicsSceneMouseEvent* event = static_cast<QGraphicsSceneMouseEvent*>(ev);
           if (event->button() == Qt::RightButton)
             this->endZoomRect(event->pos(), event->scenePos());
         }
       else if (ev->type() == QEvent::GraphicsSceneWheel)
         {
           QGraphicsSceneWheelEvent* event = static_cast<QGraphicsSceneWheelEvent*>(ev);
           int d = event->delta();
           m_auto_resize_window = false;
           m_view_center = event->scenePos();
           this->zoomBy(1.0 + d / 120 * 0.1);
           return true; // no more processing
         }

      return QObject::eventFilter(obj, ev);
    }

    /***********************************/
    /***  IMageViewer Implementation ***/
    /***********************************/

    void
    ImageViewer::_init(QtImageBase* ima)
    {
      m_ima = ima;
      m_pixmap.setPixmap(QPixmap::fromImage(m_ima->getQImage()));

      m_scene = new QGraphicsScene(this);
      QObject::connect(m_scene, SIGNAL(sceneRectChanged(const QRectF&)),
		       this, SLOT(onzoom(const QRectF&)));

      m_scene->addItem(&m_pixmap);
      m_scene->setFocus();
      m_scene->setFocusItem(&m_pixmap);

      this->setScene(m_scene);
      this->setMaximumSize(1000,1000);
      this->setHorizontalScrollBarPolicy( Qt::ScrollBarAlwaysOff );
      this->setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOff );

      m_scene->update();
      m_scene->installEventFilter(this);

      m_ev_hdler = new ImageViewerEventHandler(this);
      m_scene->installEventFilter(m_ev_hdler);
    }


    ImageViewer::~ImageViewer()
    {
      delete m_ima;
    }

    void
    ImageViewer::onzoom(const QRectF& rect)
    {
      int w = rect.right();
      int h = rect.bottom();

      this->resize(w+10, h+10);
    }

    void
    ImageViewer::reset()
    {
      m_ima->reset();
    }

    image2d<rgb8>&
    ImageViewer::getView()
    {
      return m_ima->getView();
    }

    const image2d<rgb8>&
    ImageViewer::getView() const
    {
      return m_ima->getView();
    }

    QGraphicsScene*
    ImageViewer::getScene()
    {
      return m_scene;
    }

    const QGraphicsScene*
    ImageViewer::getScene() const
    {
      return m_scene;
    }

    QGraphicsPixmapItem*
    ImageViewer::getPixmap()
    {
      return &m_pixmap;
    }

    const QGraphicsPixmapItem*
    ImageViewer::getPixmap() const
    {
      return &m_pixmap;
    }


    void
    ImageViewer::update()
    {
      std::cout << "Image updated." << std::endl;
      m_ima->update();
      m_pixmap.setPixmap(QPixmap::fromImage(m_ima->getQImage()));
    }

    void
    ImageViewer::notify()
    {
      QPixmap pxmap = m_pixmap.pixmap();
      QPainter painter(&m_ima->getQImage());
      painter.drawPixmap(pxmap.rect(), pxmap);
    }


    bool
    ImageViewer::eventFilter(QObject* obj, QEvent* ev)
    {
      if (obj == m_scene)
	{
	  if (ev->type() == QEvent::GraphicsSceneMouseMove)
	    {
	      QGraphicsSceneMouseEvent* event = static_cast<QGraphicsSceneMouseEvent*>(ev);
	      QPointF p_ = m_pixmap.mapFromScene(event->scenePos());
	      point2d p = { (short)p_.y(), (short)p_.x() };
	      if (m_ima->getView().domain().has(p)) {
		//std::cout << p << std::endl;
		emit pointHover(p);
		return false;
	      }
	    }
	  else if (ev->type() == QEvent::GraphicsSceneMousePress)
	    {
	      QGraphicsSceneMouseEvent* event = static_cast<QGraphicsSceneMouseEvent*>(ev);

              QPointF p__ = event->scenePos();
	      QPointF p_ = m_pixmap.mapFromScene(event->scenePos());
	      point2d p = { (short)p_.y(), (short)p_.x() };
	      if (m_ima->getView().domain().has(p)) {
		std::cout << "Press: " << p << "," << (point2df{p__.x(), p__.y()})
                          << std::endl;
		emit pointSelected(p);
		return false;
	      }
	    }
	}
      return QGraphicsView::eventFilter(obj, ev);
    }

    void
    ImageViewer::wheelEvent(QWheelEvent* ev)
    {
      ev->accept();
      QGraphicsView::wheelEvent(ev);
      // Disable scrolling.
    }
  }

}

