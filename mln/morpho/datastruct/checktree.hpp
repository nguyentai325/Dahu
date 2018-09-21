#ifndef MLN_MORPHO_DATASTRUCT_CHECKTREE_HPP
# define MLN_MORPHO_DATASTRUCT_CHECKTREE_HPP

# include <mln/morpho/datastruct/component_tree.hpp>

namespace mln
{
  namespace morpho
  {
    namespace internal
    {

      template <class Tree>
      void
      checktree(const Tree& tree)
      {
        typename Tree::node_type prec = tree.nend();

        mln_foreach(auto node, tree.nodes())
          {
            mln_assertion(node.get_parent_id() == prec.id() or
                          prec.get_next_sibling_id() == node.id());

            mln_assertion(node.next_node().get_parent_id() == node.id() or
                          node.get_next_node_id() == node.get_next_sibling_id());

            mln_assertion(node.next_sibling().get_parent_id() == node.get_parent_id() or
                          node.get_next_sibling_id() == node.parent().get_next_sibling_id());

            prec = node;
          }

      }
    }
  }
}

#endif // ! MLN_MORPHO_DATASTRUCT_CHECKTREE_HPP
