#include <mln/morpho/tos/tos.hpp>
#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/colors.hpp>
#include <boost/format.hpp>
#include "addborder.hpp"
#include "topology.hpp"

void usage(int argc, char** argv)
{
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " input.(ppm|png...) out_without_extension grain1 [grain2 [grain3]]" << std::endl
	      <<
      "Decompose the image in R, G, B channels, compute the Tos on each of them"
      "perform a grain filter indepedentally, and merge back the results to rgb space."
      "If a single lambda is given it outputs to 'out.tiff'."
      "Otherwise, it outputs filtered images as 'out-grain1.tiff', 'out-grain2.tiff'..."
	      << std::endl;
    std::abort();
  }
}

namespace mln
{

  /**
  * \brief compute the grain filter using ToS
  * \param ima Input image. \p ima must suitable for ToS computation, i.e having a constant border and pixel x2.
  * \pre the ToS must have been computed with a total Order (not a pre-order), if not the value of a node of equivalents pixels
  * is undefined
  * \return The filtered image.
  */
  template <typename V, typename T>
  image2d<V>
  grainfilter(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S,  unsigned grain)
  {

    // 1. compute area attribute
    image2d<unsigned> area;
    resize(area, K).init(0);

    area[S[0]] += K1::is_face_2(K.point_at_index(S[0]));
    for (int i = S.size() - 1; i > 0; --i) // Strict >0 non root
      {
	unsigned p = S[i];
	unsigned q = parent[p];
	if (K1::is_face_2(K.point_at_index(p)))
	  area[p] += 1;
	area[q] += area[p];
      }

    //2. filter
    image2d<V> out;
    resize(out, K);
    out[S[0]] = ima[S[0]];
    for (unsigned i : S)
      {
	point2d p = K.point_at_index(i);
	point2d pp = p / 2;

	if (area[i] > grain and K1::is_face_2(p))
	  out[i] = ima(pp);
	else
	  out[i] = out[parent[i]];
      }

    // 3: retrieve 2-faces
    image2d<V> output;
    box2d d = out.domain();
    resize(output, ima);
    copy(out | sbox2d(d.pmin, d.pmax, point2d{2,2}), output);

    return output;
  }

}



int main(int argc, char** argv)
{
  using namespace mln;

  usage(argc, argv);

  image2d<rgb8> ima;
  io::imread(argv[1], ima);

  typedef UInt<9> V;
  typedef image2d<V> I;
  I r = transform(ima, [](rgb8 v) -> V { return v[0] * 2; });
  I g = transform(ima, [](rgb8 v) -> V { return v[1] * 2; });
  I b = transform(ima, [](rgb8 v) -> V { return v[2] * 2; });

  I rr = addborder(r);
  I gg = addborder(g);
  I bb = addborder(b);

  image2d<V> Kr, Kg, Kb;
  image2d<unsigned> parentr, parentg, parentb;
  std::vector<unsigned> Sr, Sg, Sb;

  std::tie(Kr, parentr, Sr) = morpho::ToS(rr, c4);
  std::tie(Kg, parentg, Sg) = morpho::ToS(gg, c4);
  std::tie(Kb, parentb, Sb) = morpho::ToS(bb, c4);


  for (int i = 3; i < argc; ++i) {
    unsigned lambda = std::atoi(argv[i]);
    image2d<V> outr = grainfilter(rr, Kr, parentr, Sr, lambda);
    image2d<V> outg = grainfilter(gg, Kg, parentg, Sg, lambda);
    image2d<V> outb = grainfilter(bb, Kb, parentb, Sb, lambda);


    image2d<rgb8> res = transform(imzip(outr, outg, outb), [](const std::tuple<V,V,V>& x) {
	return rgb8{ (uint8) (std::get<0>(x) / 2), (uint8) (std::get<1>(x) / 2), (uint8) (std::get<2>(x) / 2) }; });

    std::string filename = boost::str(boost::format("%s-%05i.tiff") % argv[2] % lambda);
    io::imsave(res, filename.c_str());

  }
}




