#include <QApplication>

#include <mln/qt/imageviewer.hpp>
#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>

#include <QtGui>

int main(int argc, char** argv)
{
  using namespace mln;

  QApplication a(argc, argv);

  image2d<uint8> ima;
  io::imread(argv[1], ima);

  qt::ImageViewer main(ima);
  main.show();

  // QImage image("../../../img/small.pgm");
  // QGraphicsPixmapItem item(QPixmap::fromImage(image));
  // QGraphicsScene* m_scene = new QGraphicsScene;
  // m_scene->addItem(&item);

  // QGraphicsView* view = new QGraphicsView(m_scene);
  // QMainWindow win; 
  // win.setCentralWidget(view);
  // win.show();

  a.exec();
}
