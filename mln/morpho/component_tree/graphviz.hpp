#ifndef MLN_MORPHO_COMPONENT_TREE_GRAPHVIZ_HPP
# define MLN_MORPHO_COMPONENT_TREE_GRAPHVIZ_HPP

# include <mln/morpho/component_tree/component_tree.hpp>
# include <iosfwd>

namespace mln
{

  namespace morpho
  {

    template <class P, class Amap>
    void
    write_graphviz(std::ostream& os,
		   const component_tree<P, Amap>& tree);


    template <class P, class Amap>
    void
    write_graphviz(std::ostream& os,
		   const component_tree<P, Amap>& tree);


    /***********************************/
    /**   Implementation	      **/
    /***********************************/


    template <class P, class Amap>
    void
    write_graphviz(std::ostream& os,
		   const component_tree<P, Amap>& tree)
    {
      os << "digraph {" << std::endl;

      mln_foreach(auto x, tree.nodes())
	//os << "\t" << x.id() << " [label=\"" << property[x] << "\"];" << std::endl;
	{
	//std::cout << x.id() << std::endl;
	os << "\t" << x.id() << " [label=\"" << x.id() << "\"];" << std::endl;
	}

      //mln_foreach(auto x, tree.nodes_without_root())
      mln_foreach(auto x, tree.nodes())

	os << "\t" << x.id() << " -> " << x.get_parent_id() << ";" << std::endl;


      os << "}" << std::endl;
    }

  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_GRAPHVIZ_HPP
