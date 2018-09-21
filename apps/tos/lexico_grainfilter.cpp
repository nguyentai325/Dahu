#include <mln/morpho/tos/tos.hpp>
#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/lsh.hpp>
#include <boost/format.hpp>
#include "addborder.hpp"
#include "topology.hpp"

void usage(int argc, char** argv)
{
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << "[--lsh] input.(ppm|png...) out_without_extension grain1 [grain2 [grain3]]" << std::endl
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
  grainfilter(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S, unsigned grain)
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

    auto equiv = [](const T& x, const T& y) { return !(x < y) and !(y < x); };
    for (unsigned i : S)
      {
	point2d p = K.point_at_index(i);
	point2d pp = p / 2;
	if (equiv(K[parent[i]], K[i]))
	  area[i] = area[parent[i]];

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

  template <typename Vec>
  bool
  lex_compare(const Vec& x, const Vec& y)
  {
    for (unsigned i = 0; i < Vec::ndim; ++i)
      {
	if (x[i] < y[i])
	  return true;
	else if (x[i] > y[i])
	  return false;
      }
    return true;
  }

}


int main(int argc, char** argv)
{
  using namespace mln;

  usage(argc, argv);

  bool lsh = false;
  int pos = 1;
  if (argv[1] == std::string("--lsh")) {
    lsh = true;
    ++pos;
  }

  image2d<rgb8> ima_;
  io::imread(argv[pos], ima_);

  image2d<rgb8> ima = addborder(ima_, [] (const rgb8& x, const rgb8& y) { return x.as_vec() < y.as_vec(); });

  typedef vec<UInt<9>, 3> V; // vec are comparable by default with a lexocgraphical order
  typedef image2d<V> I;
  I f;
  if (not lsh) {
    f = transform(ima, [](rgb8 v) -> V { V x(v); x[2] *= 2; return x; }); // Note we on need the third dim to x2
  } else {
    f = transform(ima, [](rgb8 v) -> V { V x(rgb2lsh(v)); x[2] *= 2; return x; }); // Note we on need the third dim to x2
  }



  image2d<V> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  std::tie(K, parent, S) = morpho::ToS(f, c4);



  for (int i = pos+2; i < argc; ++i) {
    unsigned lambda = std::atoi(argv[i]);
    image2d<V> out = grainfilter(f, K, parent, S, lambda);

    image2d<rgb8> res;
    if (not lsh) {
      res = transform(out, [] (V v) -> rgb8 { v[2] /= 2; return (rgb8)v; });
    } else {
      res = transform(out, [] (V v) -> rgb8 { v[2] /= 2; return lsh2rgb((lsh8)v); });
    }

    std::string filename = boost::str(boost::format("%s-%05i.tiff") % argv[pos+1] % lambda);
    io::imsave(res, filename.c_str());

  }
}

