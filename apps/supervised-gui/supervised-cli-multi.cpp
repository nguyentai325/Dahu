#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/colors/literal.hpp>

#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/accu/accumulators/mean.hpp>
#include <mln/accu/accumulators/max.hpp>

#include <mln/morpho/algebraic_filter.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/croutines.hpp>
#include "constants.hpp"
#include "compute_distance.hpp"
#include <Eigen/Dense>

using namespace mln;

property_map<tree_t, uint8>
labelize_tree(const tree_t& tree,
              const image2d<uint8>& markers,
              int nlabels)
{
  typedef Eigen::VectorXi MyVec;

  property_map<tree_t, uint8> colormap(tree, NONE);
  property_map<tree_t, MyVec > count(tree, MyVec::Zero(nlabels+1) );
  mln_foreach(auto px, markers.pixels())
    {
      tree_t::node_type node = tree.get_node_at(px.index());
      if (px.val() != 0)
        count[node][px.val()]++;
    }

  mln_foreach(auto x, tree.nodes()) {
    MyVec::Index maxpos;
    int maxv = count[x].maxCoeff(&maxpos);
    if (maxv != 0)
      colormap[x] = maxpos;
  }

  return colormap;
}


image2d<uint8>
segmentation(const tree_t& tree,
             const property_map<tree_t,rgb<int> >& vmap,
             const image2d<uint8>& markers_,
             int nlabel,
             float reject = 0.5)
{


  image2d<uint8> markers = Kadjust_to(markers_, tree._get_data()->m_pmap.domain(), "zero");

   // 1. Labelize the tree
  property_map<tree_t, uint8> colormap = labelize_tree(tree, markers, nlabel);

  // 2. Compute distances
  property_map<tree_t, float> dmap[nlabel+1];

  for (int i = 1; i <= nlabel; ++i)
    dmap[i] = compute_distance(tree, colormap, vmap, i);


  // 3. Classify
  property_map<tree_t, uint8> cmap(tree, 0);
  mln_foreach(auto x, tree.nodes())
    {
      int imin = 1;
      for (int i = 1; i <= nlabel; ++i)
        if (dmap[i][x] < dmap[imin][x])
          imin = i;
      cmap[x] = imin;
    }

  image2d<uint8> out;
  resize(out, tree._get_data()->m_pmap);

  morpho::reconstruction(tree, cmap, out);
  auto output = Kadjust_to(out, markers_.domain());

  return output;
}

image2d<rgb8>
recons(const image2d<uint8>& lbl,
       const image2d<rgb8>& ori,
       int nlabel)
{
  // 4. Reconstruction
  accu::accumulators::mean<rgb8> m[nlabel+1];
  mln_pixter(pxm, pxv, lbl, ori);
  mln_forall(pxm, pxv)
    m[pxm->val()].take(pxv->val());
  m[0].init();

  rgb8 vmean[nlabel+1];
  for (int i = 0; i <= nlabel; ++i)
    vmean[i] = m[i].to_result();

  image2d<rgb8> rec = transform(lbl, [&vmean](uint8 x){ return vmean[x]; });
  return rec;
}


int main(int argc, char** argv)
{
  if (argc < 6) {
    std::cout << "Usage: " << argv[0] << " tree input.ppm markers.pgm grain output.pgm [output.ppm]\n"
      "Perform the classification on the tree and outputs a label map"
      "and a map with the average color.\n";
    std::exit(1);
  }

  tree_t tree;
  {
    std::ifstream fs(argv[1]);
    morpho::load(fs, tree);
  }

  image2d<rgb8> ima;
  io::imread(argv[2], ima);

  image2d<uint8> markers;
  io::imread(argv[3], markers);

  int grain = std::atoi(argv[4]);
  if (grain != 0) {
    grain_filter_inplace(tree, grain);
  }

  image2d<rgb8> F = Kadjust_to(ima, tree._get_data()->m_pmap.domain());
  property_map<tree_t, rgb<int> > vmap = morpho::vaccumulate_proper(tree, F, accu::features::mean<>());

  image2d<uint8> res;
  image2d<rgb8> seg;
  int nlabel = accumulate(markers, accu::features::max<> ());
  res = segmentation(tree, vmap, markers, nlabel);
  for (int i = 1; i <= nlabel; ++i) {
    auto clo = morpho::area_closing(res == i, c4, grain);
    fill(res | clo, i);
  }

  seg = recons(res, ima, nlabel);

  io::imsave(res, argv[5]);
  if (argc >= 7) {
    io::imsave(seg, argv[6]);
  }
}
