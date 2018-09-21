#ifndef BRUSH_HPP
# define BRUSH_HPP

# include <mln/qt/imageviewer.hpp>

class MyBrush : public QObject
{
  Q_OBJECT;
public:
  typedef std::function<mln::image2d<mln::rgb8>(const mln::image2d<mln::rgb8>&, float, unsigned)> callback_fun_t;

  MyBrush(mln::qt::ImageViewer* viewer,
          callback_fun_t callback);
  bool eventFilter(QObject* obj, QEvent* ev);

public slots:
  void setColor1();
  void setColor2();
  void setRadius(int k);
  void run();
  void reload();
  void revert();
  void set_reject_value(const QString& v);
  void set_ambiguity_policy(int x);

private:
  mln::qt::ImageViewer*  m_viewer;
  callback_fun_t         m_callback;
  QGraphicsScene*        m_scene;
  QGraphicsPixmapItem*   m_pixmap;
  QPixmap                m_ima;
  mln::image2d<mln::rgb8>  m_ori;
  bool                   m_active;
  QColor                 m_color;
  int                    m_radius;
  float                  m_reject_value;
  int                    m_policy;
};

#endif //!BRUSH_HPP
