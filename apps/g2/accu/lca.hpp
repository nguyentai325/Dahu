#ifndef ACCU_LCA_HPP
# define ACCU_LCA_HPP

# include <mln/accu/accumulator.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/datastruct/attribute_map.hpp>
# include <apps/tos/croutines.hpp>

namespace mln
{

  namespace accu
  {

    // Not that the argument of the accumulator is a point (or an index), not a node.
    template <class P, class Amap>
    struct least_common_ancestor : Accumulator< least_common_ancestor<P, Amap> >
    {
    private:
      typedef morpho::component_tree<P, Amap> tree_t;
      typedef property_map<tree_t, unsigned>  amap_t;

    public:
      typedef P					argument_type;
      typedef typename tree_t::node_type	result_type;

      least_common_ancestor(const tree_t& tree,
			    const amap_t& depth)
	: m_tree (tree),
	  m_depth (depth),
          m_current (m_tree.nend())
      {
      }

      void init()
      {
	m_current = m_tree.nend();
      }

      void take(const argument_type& x)
      {
	result_type xnode = m_tree.get_node_at(x);
	m_current = lca(m_tree, m_depth, m_current, xnode);
      }

      void take(const least_common_ancestor& other)
      {
	m_current = lca(m_tree, m_depth, m_current, other.m_current);
      }

      result_type to_result() const
      {
	return m_current;
      }

    private:
      const tree_t& m_tree;
      const amap_t& m_depth;

      result_type m_current;
    };

  }

}

#endif // ! ACCU_LCA_HPP
