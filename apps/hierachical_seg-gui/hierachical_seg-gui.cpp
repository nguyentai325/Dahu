
#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/literal.hpp>

#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/maxtree/maxtree.hpp>
#include <mln/data/stretch.hpp>

#include <QApplication>
#include <QtGui>
#include <mln/qt/imageviewer.hpp>
#include <apps/tos/topology.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include "brush.hpp"


using namespace mln;

typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;

///
///
/// \param tree The component tree of the saliency map
/// \param f The original image to filter out
/// \param markers__ The image with markers
/// \param[out] out The image to store the result in.
void
segmentation(const tree_t& tree,
             const image2d<rgb8>& f,
             const image2d<rgb8>& markers__,
             image2d<rgb8>& out)
{
  enum colortag { BLANC = 0, ROUGE = 1, NOIR = 2 };

  auto markers = imtransform(markers__,
                             [](const rgb8& v) -> uint8 {
                               if (v == colors::literal::red) return NOIR;
                               else if (v == colors::literal::blue) return ROUGE;
                               else return BLANC;
                             });

  property_map<tree_t, uint8> tags(tree, BLANC);
  mln_foreach(auto px, markers.pixels())
    {
      tree_t::node_type x = tree.get_node_at(px.index());
      colortag v = (colortag)px.val();

      if (K1::is_face_2(px.point()) and v != BLANC)
        {
          if (tags[x] == BLANC)
            tags[x] = v;
          else if (tags[x] != v)
            tags[x] = NOIR; // Ambuiguity => set background
        }
    }

  tags[tree.get_root()] = NOIR; // Root is background


  {
    image2d<uint8> markers = imchvalue<uint8>(markers__);
    morpho::reconstruction(tree, tags, markers);
    io::imsave(markers, "/tmp/markers.tiff");
  }


  // Propagate up
  mln_reverse_foreach(auto x, tree.nodes()) {
    if (tags[x] != BLANC) {
      auto q = x.parent();
      if (tags[q] == BLANC) // Ok propagate color
        tags[q] = tags[x];
      else if (tags[q] != tags[x]) // Either bg thus both red/black => thus black
        tags[q] = NOIR;
    }
  }

  // Propagate down
  mln_foreach(auto x, tree.nodes_without_root()) {
    auto q = x.parent();
    mln_assertion(tags[q] != BLANC);
    if (tags[x] == BLANC)
      tags[x] = tags[q];
  }

  {
    image2d<uint8> markers = imchvalue<uint8>(markers__);
    morpho::reconstruction(tree, tags, markers);
    io::imsave(markers, "/tmp/markers2.tiff");
  }

  //
  {
    mln_pixter(pxin, pxout, f, out);
    mln_forall(pxin, pxout)
      {
        auto x = tree.get_node_at(pxin->index());
        if (tags[x] == NOIR)
          pxout->val() = 192 + pxin->val() / 4;
        else
          pxout->val() = pxin->val();
      }
  }


}

int main(int argc, char** argv)
{
  QApplication a(argc, argv);

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " saliency(float) original.ppm" << std::endl;
    std::exit(1);
  }

  image2d<float> sal;
  io::imread(argv[1], sal);

  tree_t tree = morpho::mintree_indexes(sal, c4);


  image2d<rgb8> ima;
  io::imread(argv[2], ima);

  auto sal2 = data::stretch<uint8>(sal);
  ima = Kadjust_to(ima, sal.domain());

  mln_pixter(px1, px2, ima, sal2);
  mln_forall(px1, px2)
    if (px2->val())
      px1->val() = px2->val();


  QMainWindow win;
  qt::ImageViewer* viewer = new qt::ImageViewer(ima);

  auto callback = std::bind(segmentation, tree, ima,
                            std::placeholders::_1,
                            std::placeholders::_2);

  QGraphicsScene* scene = viewer->getScene();
  MyBrush brush(viewer, callback);
  scene->installEventFilter(&brush);

  QToolBar* toolbar = new QToolBar(&win);
  QAction* action1 = toolbar->addAction("Bg");
  QAction* action2 = toolbar->addAction("Fg");
  QAction* action3 = toolbar->addAction("Run");
  QAction* action4 = toolbar->addAction("Revert");
  QAction* action5 = toolbar->addAction("Reload");
  QObject::connect(action1, SIGNAL(triggered()),
                   &brush, SLOT(setColor1()));
  QObject::connect(action2, SIGNAL(triggered()),
                   &brush, SLOT(setColor2()));
  QObject::connect(action3, SIGNAL(triggered()),
                   &brush, SLOT(run()));
  QObject::connect(action4, SIGNAL(triggered()),
                   &brush, SLOT(revert()));
  QObject::connect(action5, SIGNAL(triggered()),
                   &brush, SLOT(reload()));


  QActionGroup agroup(&win);
  agroup.addAction(action1);
  agroup.addAction(action2);
  action1->setCheckable(true);
  action2->setCheckable(true);


  win.setCentralWidget(viewer);
  win.addToolBar(toolbar);
  win.show();

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
