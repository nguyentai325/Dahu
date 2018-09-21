#include "profiles.hpp"
#include <apps/g2/compute_ctos.hpp>
#include <apps/tos/croutines.hpp>
#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 3)
    {
      std::cerr << "Usage: " << argv[0] << " input.tiff stem\n";
      std::exit(1);
    }

  image2d<rgb16> f;
  io::imread(argv[1], f);

  tree_t tree = compute_ctos(f);
  grain_filter_inplace(tree, 10);

  image2d<rgb16> F;
  F = Kadjust_to(f, tree._get_data()->m_pmap.domain());

  auto vmap = morpho::pixaccumulate_proper(tree, F, mymean<rgb16>());
  {
    mln_foreach(auto x, tree.nodes_without_root()) {
      if (vmap[x] == colors::literal::black)
        vmap[x] = vmap[x.parent()];
    }
  }


  std::vector<float> lambdas;
  std::vector< image2d<rgb<int>> > recs;
  std::vector< image2d<float> > profiles = profile(tree, F, vmap, lambdas, recs);

  std::string stem = argv[2];
  int i = 0;
  for (const image2d<float>& x : profiles) {
    auto res = Kadjust_to(x, f.domain());
    io::imsave(res, (boost::format("%s-%.02f.tiff") % stem % lambdas[i]).str());
    i++;
  }
}
