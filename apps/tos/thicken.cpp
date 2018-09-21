#include <mln/core/grays.hpp>
#include <mln/core/image/image2d.hpp>
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
#include "thicken.hpp"
#include "addborder.hpp"



void usage(int argc, char** argv)
{
  if (argc < 5 or (argv[1] != std::string("mintree") && argv[1] != std::string("tos")) or
      (argv[2] != std::string("min") && argv[2] != std::string("max")))
    {
      std::cerr << "Usage: " << argv[0] << "(mintree|tos) (min|max) ima.(ppm|png|tiff...) out.tiff" << std::endl
		<< "Compute the ToS marginally, compute area attribute and compute another ToS on it."
		<< "Output the thickening" << std::endl;
      abort();
    }
}




namespace mln
{

  template <typename V, typename T, class FilterFun>
  image2d<V>
  setmean_on_nodes(const image2d<V>& ima, const image2d<T>& K,
		   const image2d<unsigned>& parent, const std::vector<unsigned>& S, FilterFun pred)
  {
    typedef rgb<unsigned> SumType;
    image2d<SumType>  sum;
    image2d<unsigned> count;
    resize(count, K).init(0);
    resize(sum, K).init(SumType ());

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
    resize(out, ima);
    unsigned p = S[0];
    out[p] = sum[p] / count[p];
    for (unsigned p: S)
      {
	unsigned q = parent[p];
	if (K[p] == K[q])
	  out[p] = out[q];
	else
	  {
	    while (count[p] == 0)
	      p = parent[p];
	    out[p] = sum[p] / count[p];
	  }
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


  auto r_area = morpho::area_compute(rK, rparent, rS);//, K1::is_face_2);
  auto g_area = morpho::area_compute(gK, gparent, gS);//, K1::is_face_2);
  auto b_area = morpho::area_compute(bK, bparent, bS);//, K1::is_face_2);

  image2d<unsigned> area;
  if (std::string(argv[2]) == "max")
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
	return std::max(std::get<0>(x), std::max(std::get<1>(x), std::get<2>(x))); });
  else
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
	return std::min(std::get<0>(x), std::min(std::get<1>(x), std::get<2>(x))); });

  //++io::imprint(area);
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


  auto w = thicken_tdn(K, parent, S, c8);
  image2d<rgb8> out;
  if (use_tos)
    out = setmean_on_nodes(tmp, w, parent, S, K2::is_face_2);
  else
    out = setmean_on_nodes(tmp, w, parent, S, K1::is_face_2);

  io::imsave(out, argv[4]);
}
