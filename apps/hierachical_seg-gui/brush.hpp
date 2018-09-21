#ifndef BRUSH_HPP
# define BRUSH_HPP

# include <mln/qt/imageviewer.hpp>

class MyBrush : public QObject
{
  Q_OBJECT;
public:
  MyBrush(mln::qt::ImageViewer* viewer,
          std::function<void(const mln::image2d<mln::rgb8>&, mln::image2d<mln::rgb8>&)> callback);
  bool eventFilter(QObject* obj, QEvent* ev);

public slots:
  void setColor1();
  void setColor2();
  void setRadius(int k);
  void run();
  void reload();
  void revert();

private:
  mln::qt::ImageViewer*  m_viewer;
  std::function<void(const mln::image2d<mln::rgb8>&,mln::image2d<mln::rgb8>&)> m_callback;
  QGraphicsScene*        m_scene;
  QGraphicsPixmapItem*   m_pixmap;
  QPixmap                m_ima;
  mln::image2d<mln::rgb8>  m_ori;
  bool                   m_active;
  QColor                 m_color;
  int                    m_radius;
};

#endif //!BRUSH_HPP
