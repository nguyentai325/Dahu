#ifndef MLN_MORPHO_COMPONENT_TREE_COMPUTE_DEPTH_HPP
# define MLN_MORPHO_COMPONENT_TREE_COMPUTE_DEPTH_HPP

# include <mln/core/trace.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/datastruct/attribute_map.hpp>

namespace mln
{

  namespace morpho
  {

    template <class P, class Amap>
    property_map<component_tree<P, Amap>, unsigned>
    compute_depth(const component_tree<P, Amap>& tree);


    /****************************/
    /**   Imlpementation       **/
    /****************************/

    template <class P, class Amap>
    property_map<component_tree<P, Amap>, unsigned>
    compute_depth(const component_tree<P, Amap>& tree)
    {
      mln_entering("mln::morpho::component_tree::compute_depth");

      property_map<component_tree<P, Amap>, unsigned> depth(tree);

      depth[tree.get_root_id()] = 0;
      mln_foreach(auto x, tree.nodes_without_root())
	depth[x.id()] = depth[x.get_parent_id()] + 1;

      mln_exiting();
      return depth;
    }

  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_COMPUTE_DEPTH_HPP
