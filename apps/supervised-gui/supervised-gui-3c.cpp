#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/colors/literal.hpp>

#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/accu/accumulators/mean.hpp>

#include <QApplication>
#include <QtGui>
#include <mln/qt/imageviewer.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include "brush.hpp"
#include "constants.hpp"
#include "compute_distance.hpp"

using namespace mln;

// Global
image2d<rgb8> original;

property_map<tree_t, uint8>
labelize_tree(const tree_t& tree,
              const image2d<uint8>& markers,
              unsigned labelization_policy)
{
  property_map<tree_t, uint8> colormap(tree, NONE);

  if ((labelization_policy & AMBIGUITY_MASK) != AMBIGUITY_AS_MAJORITY)
    {
      int ambiguity_policy = labelization_policy & AMBIGUITY_MASK;
      uint8 dummy_label;
      switch (ambiguity_policy) {
        case AMBIGUITY_AS_FOREGROUND: dummy_label = FOREGROUND; break;
        case AMBIGUITY_AS_BACKGROUND: dummy_label = BACKGROUND; break;
        default: dummy_label = DUMMY;
      }

      mln_foreach(auto px, markers.pixels())
        {
          tree_t::node_type node = tree.get_node_at(px.index());
          if (px.val() != NONE) {
            if (colormap[node] == NONE)
              colormap[node] = px.val();
            else if (colormap[node] != px.val()) // We have an ambiguity
              colormap[node] = dummy_label;
          }
        }
    }
  else
    {
      property_map<tree_t, std::pair<unsigned, unsigned> > count(tree, std::make_pair(0,0));
      mln_foreach(auto px, markers.pixels())
        {
          tree_t::node_type node = tree.get_node_at(px.index());
          if (px.val() == BACKGROUND)
            ++(count[node].first);
          else if (px.val() == FOREGROUND)
            ++(count[node].second);
        }
      mln_foreach(auto x, tree.nodes())
        if (count[x].first > count[x].second)
          colormap[x] = BACKGROUND;
        else if (count[x].first < count[x].second)
          colormap[x] = FOREGROUND;
        else
          colormap[x] = NONE;
    }

  return colormap;
}


image2d<rgb8>
segmentation(const tree_t& tree,
             const property_map<tree_t,rgb<int> >& vmap,
             const image2d<rgb8>& ori,
             const image2d<rgb8>& markers__,
             float reject = 0.5,
             unsigned policy = 0)
{
  mln_entering("Running the classification");

  image2d<uint8> markers_ = transform(markers__,
                                      [](const rgb8& v) -> uint8 {
                                        if (v == colors::literal::red) return BACKGROUND;
                                        else if (v == colors::literal::blue) return FOREGROUND;
                                        else return NONE;
                                      });

  //image2d<rgb8>  ima = Kadjust_to(ima_, tree._get_data()->m_pmap.domain());
  image2d<uint8> markers = Kadjust_to(markers_, tree._get_data()->m_pmap.domain(), "zero");


  // 1. Labelize the tree
  property_map<tree_t, uint8> colormap = labelize_tree(tree, markers, policy);

  // 2. Compute distances
  property_map<tree_t, float> d_fg = compute_distance(tree, colormap, vmap, FOREGROUND);
  property_map<tree_t, float> d_bg = compute_distance(tree, colormap, vmap, BACKGROUND);
  property_map<tree_t, float> d_dummy;

  if ((policy & AMBIGUITY_MASK) == AMBIGUITY_AS_DUMMY)
    d_dummy = compute_distance(tree, colormap, vmap, DUMMY);
  else
    d_dummy = property_map<tree_t, float> (tree, value_traits<float>::max());

  // 3. Classify

  // distance map in the order BACKGROUND, FOREGROUND, DUMMY
  property_map<tree_t, rgb<float> > dmap(tree);
  mln_foreach(auto x, tree.nodes()) {
    if (d_dummy[x] != value_traits<float>::max()) {
      float s = d_fg[x] + d_bg[x] + d_dummy[x];
      dmap[x] = {d_bg[x] / s, d_fg[x] / s, d_dummy[x] / s};
    } else {
      float s = d_fg[x] + d_bg[x];
      dmap[x] = {d_bg[x] / s, d_fg[x] / s, 1};
    }
  }

  // classification map
  // if no label has the reject threshold, it gets REJECTED
  auto cmap = make_functional_property_map<tree_t::node_type> 
    ([&dmap, reject] (const tree_t::node_type& x) -> uint8 {
      float m = minimum(dmap[x]);
      if (m > reject)
        return REJECTED;
      else if (dmap[x][0] == m)
        return BACKGROUND;
      else if (dmap[x][1] == m)
        return FOREGROUND;
      else
        return DUMMY;
    });

  image2d<rgb8> out;
  resize(out, tree._get_data()->m_pmap);

  if (trace::verbose)
    {
      {
        image2d<rgb8> out = clone(original);
        fill(out | (markers_ == BACKGROUND), colors::literal::red);
        fill(out | (markers_ == FOREGROUND), colors::literal::blue);
        std::cout << "See /tmp/markers_.tiff for the image marked" << std::endl;
        io::imsave(out, "/tmp/markers_.tiff");
      }

      {
        image2d<uint8> out;
        resize(out, tree._get_data()->m_pmap);
        morpho::reconstruction(tree, colormap, out);
        std::cout << "See /tmp/markers.tiff for tree markers" << std::endl;
        io::imsave(out, "/tmp/markers.tiff");
      }

      {
        image2d< rgb<float> > out;
        resize(out, tree._get_data()->m_pmap);
        morpho::reconstruction(tree, dmap, out);
        std::cout << "See /tmp/distances.tiff for distances" << std::endl;
        io::imsave(out, "/tmp/distances.tiff");
      }

       {
         auto bmap = make_functional_property_map<tree_t::node_type>
           ([&cmap] (const tree_t::node_type& x) -> bool { return cmap[x] == FOREGROUND; });
         image2d<bool> out;
         resize(out, tree._get_data()->m_pmap);
         morpho::reconstruction(tree, bmap, out);
         out = Kadjust_to(out, markers__.domain());
         std::cout << "See /tmp/mask.pbm for mask" << std::endl;
         io::imsave(out, "/tmp/mask.pbm");
       }
    }


  {
    mln_pixter(pxin, pxout, ori, out);
    mln_forall(pxin, pxout) {
      auto x = tree.get_node_at(pxin->index());
      uint8 color = cmap[x];
      if (color == REJECTED)
        pxout->val() = rgb8{0,0,128};
      else if (color == FOREGROUND)
        pxout->val() = pxin->val();
      else if (color == BACKGROUND)
        pxout->val() = 0;
      else
        pxout->val() = pxin->val() / 2 + rgb8{128,0,0};
    }
  }

  auto output = Kadjust_to(out, markers__.domain());

  if (trace::verbose)
    {
      std::cout << "See /tmp/segmented.png for distances" << std::endl;
      io::imsave(output, "/tmp/segmented.png");
    }

  mln_exiting();
  return output;
}


int main(int argc, char** argv)
{
  QApplication a(argc, argv);

  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << "tree input.ppm\n"
              << "If TRACE is enabled, it outputs /tmp/distances.tiff, /tmp/markers.tiff, /tmp/segmented.png"
              << std::endl;
    std::exit(1);
  }

  tree_t tree;
  {
    std::ifstream fs(argv[1]);
    morpho::load(fs, tree);
  }



  io::imread(argv[2], original);

  QMainWindow win;
  qt::ImageViewer* viewer = new qt::ImageViewer(original);


  image2d<rgb8> F = Kadjust_to(original, tree._get_data()->m_pmap.domain());
  property_map<tree_t, rgb<int> > vmap = morpho::vaccumulate_proper(tree, F, accu::features::mean<>());

  std::function<image2d<rgb8>(const image2d<rgb8>&, float, unsigned)> callback =
    std::bind(segmentation, tree, vmap, F,
              std::placeholders::_1,
              std::placeholders::_2,
              std::placeholders::_3);

  QGraphicsScene* scene = viewer->getScene();
  MyBrush brush(viewer, callback);
  scene->installEventFilter(&brush);

  QToolBar* toolbar = new QToolBar(&win);
  QLineEdit* threshold_widget = new QLineEdit("0.5");
  QComboBox* ambiguity_policy_selector = new QComboBox;
  ambiguity_policy_selector->addItem("FOREGROUND", QVariant(AMBIGUITY_AS_FOREGROUND));
  ambiguity_policy_selector->addItem("BACKGROUND", QVariant(AMBIGUITY_AS_BACKGROUND));
  ambiguity_policy_selector->addItem("3rd CLASS", QVariant(AMBIGUITY_AS_DUMMY));
  ambiguity_policy_selector->addItem("MAJORITY", QVariant(AMBIGUITY_AS_MAJORITY));
  ambiguity_policy_selector->addItem("UNTAGGED", QVariant(AMBIGUITY_AS_UNTAGGED));

  QAction* action1 = toolbar->addAction("Bg");
  QAction* action2 = toolbar->addAction("Fg");
  QAction* action3 = toolbar->addAction("Run");
  QAction* action4 = toolbar->addAction("Revert");
  QAction* action5 = toolbar->addAction("Reload");
  toolbar->addWidget(ambiguity_policy_selector);
  toolbar->addWidget(threshold_widget);


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
  QObject::connect(threshold_widget, SIGNAL(textEdited(const QString&)),
                   &brush, SLOT(set_reject_value(const QString&)));
  QObject::connect(ambiguity_policy_selector, SIGNAL(currentIndexChanged(int)),
                   &brush, SLOT(set_ambiguity_policy(int)));

  QActionGroup agroup(&win);
  agroup.addAction(action1);
  agroup.addAction(action2);
  action1->setCheckable(true);
  action2->setCheckable(true);


  win.setCentralWidget(viewer);
  win.addToolBar(toolbar);
  win.show();

  a.exec();

}
