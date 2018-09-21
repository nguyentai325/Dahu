#ifndef MLN_QT_IMAGEVIEWER_HPP
# define MLN_QT_IMAGEVIEWER_HPP

# include <QGraphicsScene>
# include <QGraphicsView>
# include <QGraphicsPixmapItem>
# include <mln/qt/qtimage.hpp>

namespace mln
{
  namespace qt
  {


    class ImageViewerEventHandler; // Fwd declaration

    class ImageViewer : public QGraphicsView
    {
      Q_OBJECT;

    public:
      template <typename V>
      ImageViewer(const image2d<V>& ima, QWidget* parent = 0);

      virtual ~ImageViewer();

      void reset();
      image2d<rgb8>&		getView();
      const image2d<rgb8>&	getView() const;
      QGraphicsScene*           getScene();
      const QGraphicsScene*     getScene() const;
      QGraphicsPixmapItem*           getPixmap();
      const QGraphicsPixmapItem*     getPixmap() const;


    private:
      void _init(QtImageBase* ima);

    private:
      friend class ImageViewerEventHandler;

      QGraphicsPixmapItem                       m_pixmap;
      QGraphicsScene*                           m_scene;
      QtImageBase*                              m_ima;
      ImageViewerEventHandler*                  m_ev_hdler;

    protected:
      virtual bool eventFilter(QObject* obj, QEvent* event);
      virtual void wheelEvent(QWheelEvent* e);

    signals:
      void pointSelected(const mln::point2d& pt);
      void pointHover(const point2d& pt);

    public slots:
      /// \brief Update the on-screen pixmap image from the underlying
      /// milena's image.
      void update();

      /// \brief Write the on-screen pixmap image to the underlying milena's
      /// image (accessible through \p getView())
      void notify();

    protected slots:
      void onzoom(const QRectF& rect);
    };


    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    template <typename V>
    ImageViewer::ImageViewer(const image2d<V>& ima, QWidget* parent)
      : QGraphicsView(parent)
    {
      _init(new QtImage<V>(ima));
    }



  } // end of namespace mln::qt

} // end of namespace mln

#endif //!MLN_QT_IMAGEVIEWER_HPP
