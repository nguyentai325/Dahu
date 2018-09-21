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
#include <boost/program_options.hpp>
#include "topology.hpp"
#include <mln/morpho/maxtree_ufind_parallel.hpp>
#include <libgen.h>
#include "set_mean_on_nodes.hpp"
#include "addborder.hpp"
#include "gradient.hpp"


namespace po = boost::program_options;


po::variables_map
usage(int argc, char** argv)
{

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("close", po::value<int>(), "Perform a closing.");

  po::options_description hidden("Allowed options");
  hidden.add_options()
    ("tree", po::value<std::string>()->required(), "Type of tree (mintree|tos).")
    ("merge", po::value<std::string>()->required(), "Merge policy (min|max).")
    ("input", po::value<std::string>()->required(), "Input color file (ppm|png...)")
    ("output", po::value<std::string>()->required(), "Output stem (without extension")
    ("grain", po::value< std::vector<int> >()->composing(), "Grain filter parameters.")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::positional_options_description p;
  p.add("tree", 1);
  p.add("merge", 1);
  p.add("input", 1);
  p.add("output", 1);
  p.add("grain", -1);

  po::variables_map vm;

  try {
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
      throw;
    std::string tree = vm["tree"].as<std::string> ();
    std::string merge = vm["merge"].as<std::string> ();
    if (tree != "mintree" and tree != "tos")
      throw ;
    if (merge != "min" and merge != "max" and merge != "comp" and merge != "grad")
      throw ;
  } catch (...) {
    std::cout << "Usage: " <<  argv[0]
	      << " [options] (tos|mintree) (min|max|comp|grad) input.(tiff|ppm...) out(without extension) [lambda1 [lambda2[...]]]" << std::endl
	      << "Compute the ToS marginally, compute area attribute and compute another ToS on it." << std::endl
	      << "Output the results of filtering" << std::endl
	      << desc << std::endl;
    throw;
  }

  return vm;
}


namespace mln
{

  template <typename V>
  void
  close(image2d<V>& K, image2d<unsigned>& parent, const std::vector<unsigned>& S, unsigned lambda)
  {
    image2d<unsigned> count;
    resize(count, K).init(0);

    for (int i = S.size() - 1; i >= 0; --i)
      {
	unsigned p = S[i];
	unsigned q = parent[p];
	if (K[p] != K[q])
	  q = p;
	count[q] += 1;
      }

    for (unsigned p: S)
      {
	unsigned q = parent[p];
	if (K[p] == K[q])
	  count[p] = count[q];

	if (count[p] < lambda) {
	  K[p] = K[q];
	  if (count[q] >= lambda)
	    parent[p] = q;
	  else
	    parent[p] = parent[q];
	}
      }
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


int main(int argc, char** argv)
{
  using namespace mln;

  po::variables_map vm = usage(argc, argv);

  std::string filename = vm["input"].as<std::string> ();

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
  std::string merge = vm["merge"].as<std::string>();
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


  unsigned maxr = r_area[rS[0]];
  unsigned maxg = g_area[gS[0]];
  unsigned maxb = b_area[bS[0]];
  std::cout << maxr << " " << maxg << " " << maxb << std::endl;
  io::imsave(transform(r_area, [=](unsigned x) -> float { return (float)x / maxr; }), "red.tiff");
  io::imsave(transform(g_area, [=](unsigned x) -> float { return (float)x / maxg; }), "green.tiff");
  io::imsave(transform(b_area, [=](unsigned x) -> float { return (float)x / maxb; }), "blue.tiff");


  //++io::imprint(area);
  image2d<unsigned> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  bool use_tos = vm["tree"].as<std::string> () == "tos";

  if (use_tos)
    std::tie(K, parent, S) = morpho::ToS(area, c4);
  else {
    K = area;
    std::tie(parent, S) = morpho::impl::serial::maxtree_ufind(area, c8, std::greater<unsigned> ());
  }

  auto ima2 = addborder(ima, lexicographicalorder_less<rgb8>() ); // add border with median w.r.t < lexico
  image2d<rgb8> tmp;
  resize(tmp, parent).init(rgb8{0,0,255});

  point2d strides = use_tos ? point2d{4,4} : point2d{2,2};
  copy(ima2, tmp | sbox2d(tmp.domain().pmin, tmp.domain().pmax, strides));


  {
    image2d<unsigned> x;
    if (use_tos)
      x = morpho::area_compute(K, parent, S, K2::is_face_2);
    else
      x = morpho::area_compute(K, parent, S, K1::is_face_2, std::greater<unsigned> ());
    io::imsave(transform(x, [=](unsigned v) -> float { return (float)v; }), "area2.tiff");
  }

  // Set mean on nodes
  // if (use_tos)
  //   tmp = set_mean_on_node(tmp, K, S, parent, K2::is_face_2);
  // else
  //   tmp = set_mean_on_node(tmp, K, S, parent, K1::is_face_2);

  if (vm.count("grain"))
    {
      std::vector<int> lambdas = vm["grain"].as< std::vector<int> >();
      std::string stem = vm["output"].as<std::string> ();
      for (int fvalue: lambdas) {
	image2d<rgb8> out;
	if (use_tos)
	  out = morpho::area_filter(tmp, K, parent, S, fvalue, K2::is_face_2);
	else
	  out = morpho::area_filter(tmp, K, parent, S, fvalue, K1::is_face_2, rgb<unsigned> (),  std::greater<unsigned> ());

	image2d<rgb8> under;
	resize(under, ima2);
	copy(out | sbox2d(out.domain().pmin, out.domain().pmax, strides), under);
	io::imsave(under, boost::str(boost::format("%s-%06i.tiff") % stem % fvalue).c_str());
      }
    }

}
