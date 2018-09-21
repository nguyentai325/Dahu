#include <mln/morpho/tos/immerse.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/morpho/filtering.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>
#include <libgen.h>
#include "addborder.hpp"
#include "topology.hpp"
#include <mln/morpho/maxtree_pqueue_parallel.hpp>
#include <mln/core/algorithm/fill.hpp>

namespace mln
{

  template <typename T>
  void
  displayTree(const image2d<T>& ima,
	      const image2d<typename image2d<T>::size_type>& parent,
	      const std::vector<typename image2d<T>::size_type>& S,
	      const point2d strides = point2d{1,1})
  {
    mln_precondition(ima.domain() == parent.domain());
    mln_precondition(ima.domain().size() == S.size());

    image2d<bool> deja_vu;
    point2d pmin = ima.domain().pmin;
    point2d pmax = ima.domain().pmax;


    io::imprint(ima);
    resize(deja_vu, ima);
    auto deja_vu_ = deja_vu | sbox2d(pmin, pmax, strides);

    // Print each node
    for (unsigned i = 0; i < S.size(); ++i)
      {
	fill(deja_vu, false);
	unsigned q = S[i];
	if (i != 0 and ima[q] == ima[parent[q]])
	  continue;

	std::cout << "==== showing node " << q << " @ " << ima[q] << std::endl;
	deja_vu[q] = true;
	for (unsigned j = i; j < S.size() ; ++j)
	  {
	    unsigned p = S[j];
	    if (deja_vu[parent[p]])
	      deja_vu[p] = true;
	  }
	Kdisplay( deja_vu_, strides );
      }
  }

}

void usage(int argc, char** argv)
{
  if (argc < 4 or (argv[1] != std::string("mintree") && argv[1] != std::string("tos")) or
      (argv[2] != std::string("min") && argv[2] != std::string("max")))
    {
      std::cerr << "Usage: " << argv[0] << "(mintree|tos) (min|max) ima.(ppm|png|tiff...) lambda1 [lambda2  [lambda3 [...]]]" << std::endl
		<< "Compute the ToS marginally, compute area attribute and compute another ToS on it."
		<< "Output the results of filtering" << std::endl;
      abort();
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

  auto r_area = morpho::area_compute(rK, rparent, rS, K1::is_face_2);
  auto g_area = morpho::area_compute(gK, gparent, gS, K1::is_face_2);
  auto b_area = morpho::area_compute(bK, bparent, bS, K1::is_face_2);

  image2d<unsigned> area;
  if (std::string(argv[2]) == "max")
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
	return std::max(std::get<0>(x), std::max(std::get<1>(x), std::get<2>(x))); });
  else
    area = transform(imzip(r_area, g_area, b_area), [](const std::tuple<unsigned, unsigned, unsigned>& x) {
	return std::min(std::get<0>(x), std::min(std::get<1>(x), std::get<2>(x))); });

  unsigned maxr = r_area[rS[0]];
  unsigned maxg = g_area[gS[0]];
  unsigned maxb = b_area[bS[0]];
  std::cout << maxr << " " << maxg << " " << maxb << std::endl;
  io::imsave(transform(r_area, [=](unsigned x) -> float { return (float)x / maxr; }), "red.tiff");
  io::imsave(transform(g_area, [=](unsigned x) -> float { return (float)x / maxg; }), "green.tiff");
  io::imsave(transform(b_area, [=](unsigned x) -> float { return (float)x / maxb; }), "blue.tiff");
  io::imsave(transform(area, [=](unsigned x) -> float { return (float)x / maxr; }), "area.tiff");
  // io::imsave(r_area, "red.tiff");
  // io::imsave(g_area, "green.tiff");
  // io::imsave(b_area, "blue.tiff");


  //++io::imprint(area);
  image2d<unsigned> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  bool use_tos = argv[1] == std::string("tos");

  if (use_tos)
    std::tie(K, parent, S) = morpho::ToS(area, c4);
  else
    K = area,
    std::tie(parent, S) = morpho::impl::serial::maxtree_pqueue(area, c4, std::greater<unsigned> ());

  std::cout << "S.size(): " << S.size() << std::endl;
  auto ima2 = addborder(ima); // add border with median w.r.t < lexico
  image2d<rgb8> tmp;
  resize(tmp, parent, parent.border(), rgb8{0,0,255});

  point2d strides = use_tos ? point2d{4,4} : point2d{2,2};
  copy(ima2, tmp | sbox2d(tmp.domain().pmin, tmp.domain().pmax, strides));

  if (use_tos)
    displayTree(K, parent, S, point2d{2,2});
  else
    displayTree(K, parent, S, point2d{1,1});
}
