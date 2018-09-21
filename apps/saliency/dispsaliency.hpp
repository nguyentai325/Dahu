#ifndef DISPSALIENCY_HPP
# define DISPSALIENCY_HPP

#include <QObject>
#include <QtGui>
#include <QSlider>

#include <mln/qt/imageviewer.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/accu/accumulators/minmax.hpp>

#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/image/sub_image.hpp>

#include <mln/labeling/blobs.hpp>
#include <mln/labeling/accumulate.hpp>
#include <mln/accu/accumulators/mean.hpp>

#include <mln/io/imprint.hpp>

using namespace mln;

struct Displayer : public QObject
{
  Q_OBJECT

public:

  Displayer(const image2d<float>& f,
	    const image2d<rgb8>* ref = NULL,
	    const QString& title = "")
    : m_ori(f),
      m_ref (ref),
      m_ima((f.nrows()+1)/2, (f.ncols()+1)/2),
      m_win(m_ima)
  {
    float vmin, vmax;
    std::tie(vmin, vmax) = accumulate(f, accu::features::minmax<> ());

    m_slider.setOrientation(Qt::Horizontal);
    m_slider.setMinimum(vmin * 1000);
    m_slider.setMaximum(vmax * 1000);

    std::cout << vmin << " / " << vmax << std::endl;

    QObject::connect(&m_slider, SIGNAL(sliderReleased()),
		     this, SLOT(on_filtering()));

    this->on_filtering();

    m_win.setWindowTitle(title);
    m_win.show();
    m_slider.show();
  }

  static
  rgb8 crandom_lut(unsigned x)
  {
    auto rol = [](int x, int n) {
      return (x << n) | (x > (16-n));
    };

    unsigned k = x % 6;
    std::tuple<uint8, uint8,uint8> v {(rol(x,2) & 255),
	(rol(x,5) & 255),
	(rol(x,7) & 255) };
    uint8 r,g,b;
    switch (k) {
      case 0: std::tie(r,g,b) = v; break;
      case 1: std::tie(r,b,g) = v; break;
      case 2: std::tie(g,r,b) = v; break;
      case 3: std::tie(g,b,r) = v; break;
      case 4: std::tie(b,g,r) = v; break;
      default: std::tie(b,r,g) = v; break;
    };

    return rgb8{r,g,b};
  }

public slots:

  void on_filtering()
  {
    float v = m_slider.value() / 1000.0;
    auto bin = m_ori <= v;

    image2d<unsigned> lbl;
    unsigned nlabel;

    std::tie(lbl, nlabel) = labeling::blobs(bin, c4, 0u);
    std::cout << "Label: " << nlabel << std::endl;
    // K1 -> K0

    sbox2d dom = {{0,0}, {m_ori.nrows(), m_ori.ncols()}, {2,2}};
    if (m_ref == NULL)
      {
	auto tmp =  imtransform(lbl | dom, &Displayer::crandom_lut);
	copy(tmp, m_ima);
      }
    else
      {
	auto means = labeling::v_accumulate(lbl | dom, *m_ref, nlabel, accu::features::mean<> ());
	for (auto x: means)
	  std::cout << x << std::endl;

	auto tmp = imtransform(lbl | dom, [&means](unsigned i) { return means[i]; });
	copy(tmp, m_ima);
      }

    m_win.reset();
    m_win.update();
  }



private:
  image2d<float>	m_ori;
  const image2d<rgb8>*	m_ref;
  image2d<rgb8>		m_ima;
  qt::ImageViewer	m_win;
  QSlider		m_slider;
};

#endif // ! DISPSALIENCY_HPP
