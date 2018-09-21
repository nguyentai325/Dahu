#ifndef COMPUTE_DISTANCE_HPP
# define COMPUTE_DISTANCE_HPP

# include "constants.hpp"
# include <mln/colors/lab.hpp>

typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

/// \brief Compute the distance to a label
/// \param tree The tree used as the topological space
/// \param colormap The associative map that associates a tag to each node
/// \param vmap The associative map thata ssociated a value to each node (used for distances)
/// \param color The reference tag
inline
mln::property_map<tree_t, float>
compute_distance(const tree_t& tree,
                 const mln::property_map<tree_t, mln::uint8>& colormap,
                 const mln::property_map<tree_t, mln::rgb<int> >& vmap,
                 int color)
{
  using namespace mln;

  mln_entering("Computing distance");
  // Chamfer algorithm on the tree
  property_map<tree_t, float> distancemap(tree);

  mln_foreach(auto x, tree.nodes())
    distancemap[x] = (colormap[x] == color) ? 0 : value_traits<float>::sup();

  distancemap[tree.npos()] = value_traits<float>::sup();

  //auto mydist = [](rgb8 a, rgb8 b) -> float { return l2norm(rgb2lab(a) - rgb2lab(b)); };
  auto mydist = [](rgb8 a, rgb8 b) -> float { return l2norm(a - b); };
  // Upward
  mln_reverse_foreach(auto x, tree.nodes_without_root())
    {
      // Note the 0.1
      // The distance must stricly > 0
      float tmp = distancemap[x] + mydist(vmap[x], vmap[x.parent()]) + 0.1;
      if (tmp < distancemap[x.parent()])
        distancemap[x.parent()] = tmp;
    }


  // Downward
  mln_foreach(auto x, tree.nodes_without_root())
    {
      // The distance must stricly > 0
      float tmp = distancemap[x.parent()] + mydist(vmap[x], vmap[x.parent()]) + 0.1;
      if (tmp < distancemap[x])
        distancemap[x] = tmp;
    }


  mln_exiting();
  return distancemap;
}



#endif // ! COMPUTE_DISTANCE_HPP
