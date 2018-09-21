#ifndef SATMAXTREE_HPP
# define SATMAXTREE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>

namespace mln
{

  /// \brief Compute the saturated maxtree of an image
  std::pair<
    morpho::component_tree<unsigned, image2d<unsigned> >,
    property_map<morpho::component_tree<unsigned, image2d<unsigned> >, uint16>
    >
  satmaxtree(const image2d<uint16>& f, point2d pinf = {0,0});



  /// \brief Remove the 2F in the tree.
  morpho::component_tree<unsigned, image2d<unsigned> >
  tree_keep_2F(const morpho::component_tree<unsigned, image2d<unsigned>>& tree);

}



#endif // ! SATMAXTREE_HPP

