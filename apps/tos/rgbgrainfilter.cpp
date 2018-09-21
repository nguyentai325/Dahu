#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/morpho/filtering.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>
#include "topology.hpp"
#include <mln/morpho/maxtree_ufind_parallel.hpp>
#include <libgen.h>
#include "addborder.hpp"
#include "gradient.hpp"


void usage(int argc, char** argv)
{
  std::string type, merge;
  if (argc < 5)
    goto kill;
  type = argv[1];
  if (type != "mintree" and type != "tos")
    goto kill;
  merge = argv[2];
  if (merge != "min" and merge != "max" and merge != "comp" and merge != "grad")
    goto kill;
  return;

 kill:
  std::cerr << "Usage: " << argv[0] << "(mintree|tos) (min|max) ima.(ppm|png|tiff...) out_without_extension lambda1 [lambda2  [lambda3 [...]]]" << std::endl
	    << "Compute the ToS marginally, compute area attribute and compute another ToS on it."
	    << "Output the results of filtering" << std::endl;
  abort();
}


template < typename V, typename Compare=std::less<V> >
std::tuple<V,V,V>
minmedmax(const V& x, const V& y, const V& z, const Compare& cmp = Compare())
{
  if (cmp(x,y)) {
    if (cmp(y,z))
      return std::make_tuple(x,y,z);
    else if (cmp(x,z))
      return std::make_tuple(x,z,y);
    else
      return std::make_tuple(z,x,y);
  } else {
    if (cmp(x,z))
      return std::make_tuple(y,x,z);
    else if (cmp(z,y))
      return std::make_tuple(z,y,x);
    else
      return std::make_tuple(y,z,x);
  }
}
namespace mln
{
  template <typename V, typename T, class FilterFun>
  image2d<V>
  setmean_on_nodes(const image2d<V>& ima, const image2d<T>& K,
		   const image2d<unsigned>& parent, const std::vector<unsigned>& S, FilterFun pred)
  {
    mln_precondition(ima.domain() == K.domain());
    mln_precondition(parent.domain() == K.domain());

    typedef rgb<unsigned> SumType;
    image2d<SumType>  sum;
    image2d<unsigned> count;
    resize(count, K, K.border(), 0);
    resize(sum, K, K.border(), SumType ());

    for (int i = S.size() - 1; i >= 0; --i)
      {
	unsigned p = S[i];
	if (pred(K.point_at_index(p)))
	  {
	    unsigned q = parent[p];
	    if (K[p] != K[q])
	      q = p;
	    count[q] += 1;
	    sum[q] += ima[p];
	  }
      }

    image2d<V> out;
    resize(out, ima, ima.border(), literal::zero);
    unsigned p = S[0];
    out[p] = sum[p] / count[p];
    for (unsigned p: S)
      {
	unsigned q = parent[p];
        assert(K[q] != K[parent[q]] or q == parent[q]);

	if (K[p] == K[q] or count[p] == 0) {
	  out[p] = out[q];
        } else {
	  out[p] = sum[p] / count[p];
	}
        assert(out[p] != literal::zero);
      }
    return out;
  }



  template <typename T>
  image2d<T>
  interpolate_k1(const image2d<T>& ima)
  {
    image2d<T> out(2*ima.nrows()-1, 2*ima.ncols()-1);
    typedef point2d P;
    mln_foreach(point2d p, ima.domain())
      {
	T a = ima.at(p),
	  b = ima.at(p + P{0,1}),
	  c = ima.at(p + P{1,0}),
	  d = ima.at(p + P{1,1});

	point2d q = 2 * p;
	out.at(q) = ima.at(p);
	out.at(q + P{0,1}) = (a + b) / 2;
	out.at(q + P{1,0}) = (a + c) / 2;
	out.at(q + P{1,1}) = (a + b + c + d) / 4;
      }

    return out;

  }
}


int main(int argc, char** argv)
{
  using namespace mln;

  usage(argc, argv);

  const char* filename = argv[3];

  image2d<rgb8> ima;
  io::imread(filename, ima);

  typedef UInt<9> V;
  typedef image2d<V> I;
  I r = transform(ima, [](rgb8 v) -> V { return v[0] * 2; });
  I g = transform(ima, [](rgb8 v) -> V { return v[1] * 2; });
  I b = transform(ima, [](rgb8 v) -> V { return v[2] * 2; });
  I rr = addborder(r);
  I gg = addborder(g);
  I bb = addborder(b);


  image2d<V> rK, gK, bK;
  image2d<unsigned> rparent, gparent, bparent;
  std::vector<unsigned> rS, gS, bS;

  std::tie(rK, rparent, rS) = morpho::ToS(rr, c4);
  std::tie(gK, gparent, gS) = morpho::ToS(gg, c4);
  std::tie(bK, bparent, bS) = morpho::ToS(bb, c4);

  // io::imprint(K);
  // io::imprint(parent);

  auto r_area = morpho::area_compute(rK, rparent, rS);//, K1::is_face_2);
  auto g_area = morpho::area_compute(gK, gparent, gS);//, K1::is_face_2);
  auto b_area = morpho::area_compute(bK, bparent, bS);//, K1::is_face_2);


  image2d<unsigned> area;
  std::string merge = std::string(argv[2]);
  if (merge == "max") {
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
  	return std::max(std::get<0>(x), std::max(std::get<1>(x), std::get<2>(x))); });
    io::imsave(transform(area, [=](unsigned x) -> float { return (float)x; }), "maxarea.tiff");
  } else if (merge == "min") {
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
  	return std::min(std::get<0>(x), std::min(std::get<1>(x), std::get<2>(x))); });
    io::imsave(transform(area, [=](unsigned x) -> float { return (float)x; }), "minarea.tiff");
  } else if (merge == "comp") {
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
	unsigned min, med, max;
	std::tie(min, med, max) = minmedmax(std::get<0>(x), std::get<1>(x), std::get<2>(x));
	return (med - min) < (max - med) ? min : max;
    });
    io::imsave(transform(area, [=](unsigned x) -> float { return (float)x; }), "comparea.tiff");
  } else if (merge == "grad") {
    int size = 7;
    auto grad_r = interpolate_k1(gradient(rr, size));
    auto grad_g = interpolate_k1(gradient(gg, size));
    auto grad_b = interpolate_k1(gradient(bb, size));

    resize(area, r_area);
    mln_foreach(point2d p, area.domain())
    {
      int r = grad_r(p), g = grad_g(p), b = grad_b(p);
      int ar = r_area(p), ag = g_area(p), ab = b_area(p);

      if (r >= g) {
	area(p) = (r >= b) ? ar : ab;
      } else {
	area(p) = (g >= b) ? ag : ab;
      }
    }
    io::imsave(transform(area, [=](unsigned x) -> float { return (float)x; }), "gradarea.tiff");
  }


  io::imsave(transform(area, [] (unsigned x) -> float { return x; }), "area.tiff");

  image2d<unsigned> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  bool use_tos = argv[1] == std::string("tos");

  if (use_tos)
    std::tie(K, parent, S) = morpho::ToS(area, c4);
  else
    K = area,
    std::tie(parent, S) = morpho::impl::serial::maxtree_ufind(area, c8, std::greater<unsigned> ());

  auto ima2 = addborder(ima, lexicographicalorder_less<rgb8>() ); // add border with median w.r.t < lexico
  image2d<rgb8> tmp;
  resize(tmp, parent).init(rgb8{0,0,255});

  point2d strides = use_tos ? point2d{4,4} : point2d{2,2};
  copy(ima2, tmp | sbox2d(tmp.domain().pmin, tmp.domain().pmax, strides));


  for (int i = 5; i < argc; ++i) {
    int fvalue = std::atoi(argv[i]);
    image2d<rgb8> out, out2;
    if (use_tos) {
      out = morpho::area_filter(tmp, K, parent, S, fvalue, K2::is_face_2);
      out2 = setmean_on_nodes(out, K, parent, S, K2::is_face_2);
    } else {
      out = morpho::area_filter(tmp, K, parent, S, fvalue, K1::is_face_2, rgb<unsigned> (), std::greater<unsigned> ());
      out2 = setmean_on_nodes(out, K, parent, S, K1::is_face_2);
    }
    image2d<rgb8> under, under2;
    resize(under, ima2);
    resize(under2, ima2);
    copy(out | sbox2d(out.domain().pmin, out.domain().pmax, strides), under);
    copy(out2 | sbox2d(out2.domain().pmin, out2.domain().pmax, strides), under2);
    io::imsave(out, boost::str(boost::format("%s-%06i_k1.tiff") % argv[4] % fvalue).c_str());
    io::imsave(under, boost::str(boost::format("%s-%06i.tiff") % argv[4] % fvalue).c_str());
    io::imsave(out2, boost::str(boost::format("%s_2-%06i_k1.tiff") % argv[4] % fvalue).c_str());
    io::imsave(under2, boost::str(boost::format("%s_2-%06i.tiff") % argv[4] % fvalue).c_str());
  }

}
