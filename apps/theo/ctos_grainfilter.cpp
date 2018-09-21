#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/accu/accumulators/mean.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/croutines.hpp>
#include <apps/g2/compute_ctos.hpp>

enum e_usage_mode {
  USE_TOS, USE_MINTREE, USE_MAXTREE
};

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " [TREE] input.(ppm,jpg...) grain output.(jpg,png...)\n"
              << "  TREE must be in {tos,mintree,maxtree}\n";
    std::exit(1);
  }

  std::string cmode = argv[1];
  e_usage_mode mode;
  if (cmode == "tos")
    mode = USE_TOS;
  else if (cmode == "mintree")
    mode = USE_MINTREE;
  else if (cmode == "maxtree")
    mode = USE_MAXTREE;
  else {
    std::cerr << "TREE must be in {tos,mintree,maxtree}\n";
    std::exit(1);
  }

  image2d<rgb8> f;
  io::imread(argv[2], f);

  image2d<rgb8> ima, F;
  if (mode == USE_TOS) {
    ima = addborder_marginal(f);
    F = interpolate_k1(ima);
  } else {
    ima = f;
    F = f;
  }

  typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;
  tree_t tree;
  property_map<tree_t, size_t> areamap;

  if (mode == USE_TOS) {
    tree = compute_ctos(f);

    auto isface2 = [](const point2d& p) { return K1::is_face_2(p); };
    typedef accu::accumulators::accu_if< accu::accumulators::count<>, decltype(isface2), point2d> ACC;
    ACC accu(accu::accumulators::count<> (), isface2);
    areamap = morpho::paccumulate(tree, F, accu);
  } else if (mode == USE_MINTREE) {
    tree = compute_ctos_from_maxtrees(f, NULL, true);
    areamap = morpho::paccumulate(tree, F, accu::features::count<> ());
  } else {
    tree = compute_ctos_from_maxtrees(f, NULL, false);
    areamap = morpho::paccumulate(tree, F, accu::features::count<> ());
  }

  auto vmap = morpho::vaccumulate_proper(tree, F, accu::features::mean<>());

  int lambda = std::atoi(argv[3]);

  mln_foreach(auto x, tree.nodes())
    if (areamap[x] < lambda)
      vmap[x] = vmap[x.parent()];

  mln_foreach(auto px, F.pixels())
    {
      tree_t::node_type x = tree.get_node_at(px.index());
      if (areamap[x] < lambda)
        px.val() = vmap[x];
    }

  //grain_filter_inplace(tree, grain, false);
  // image2d<rgb8> R;
  // resize(R, F);
  // morpho::reconstruction(tree, vmap, R);

  F = Kadjust_to(F, f.domain());
  io::imsave(F, argv[4]);
}
