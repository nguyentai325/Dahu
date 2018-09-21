#include <QApplication>
#include <QtGui>

#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>

#include <mln/io/imread.hpp>
#include <mln/morpho/tos/tos.hpp>
#include <mln/morpho/filtering.hpp>

#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/addborder.hpp>
//#include <apps/tos/objdetection.hpp>
//#include <apps/attributes/MSER.hpp>
#include <mln/qt/imageviewer.hpp>
#include "qattribute.hpp"
#include "dispatcher.hpp"
#include "plotwindow.hpp"
#include "attributes/area.hpp"
#include "attributes/gray.hpp"
#include "attributes/mser.hpp"
#include "attributes/meaningfullness.hpp"

#include <apps/saliency/extinction.hpp>

int main(int argc, char** argv)
{
  using namespace mln;


  if (argc < 2)
    std::cerr << "Usage: " << argv[0] << " ima.pgm xpmin ypmin" << std::endl;


  image2d<uint8> ima;
  io::imread(argv[1], ima);

  typedef UInt<9> V;
  image2d<uint8> bima = addborder(ima);
  image2d<uint8> f = interpolate_k1(bima);
  image2d<V> ima_ = transform(bima, [](uint8 x) -> V { return x*2; });

  image2d< UInt<9> > K;
  std::vector<unsigned> S;
  image2d<unsigned> parent;

  point2d pmin;
  if (argc > 2)
    {
      pmin[0] = std::atoi(argv[2]);
      pmin[1] = std::atoi(argv[3]);
    }
  else
    {
      pmin = ima_.domain().pmin;
    }

  std::tie(K, parent, S) = morpho::ToS(ima_, c4, pmin, std::less<V> (), std::equal_to<V> ());

  QApplication a(argc, argv);

  // Create image to display
  auto Kui8 = transform(K, [](const UInt<9>& v) -> uint8 {
      return v / 2;
    });
  qt::ImageViewer w1(Kui8);
  qt::ImageViewer w2(Kui8);

  w1.show();
  w2.show();

  // Create plot window
  PlotWindow win;
  AreaAttribute<V> area(K, parent, S);
  GrayLevelAttribute<V> graylevel(K, parent, S);
  MSERAttribute<V> mser(K, parent, S);
  MeaningFullnessAttribute<uint8, V> meaningfness(bima, K, parent, S);

  win.register_attribute(&area);
  win.register_attribute(&graylevel);
  win.register_attribute(&mser);
  win.register_attribute(&meaningfness);
  win.show();


  // Set dispatcher
  QDispatcher disp(K, parent, S);
  disp.addImageWindow(&w1);
  disp.addImageWindowToFilter(&w2,
			      transform(Kui8, [] (uint8 x) -> rgb8 { return rgb8{x,x,x}; }));

  disp.setPlotWindow(&win);

  a.exec();
}
