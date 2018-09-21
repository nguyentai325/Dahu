#ifndef MLN_MORPHO_COMPONENT_TREE_ATTRIBUTE_MAP_HPP
# define MLN_MORPHO_COMPONENT_TREE_ATTRIBUTE_MAP_HPP

# include <mln/core/property_map.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>

namespace mln
{


    template <class P, class AssociativeMap, class V>
    struct property_map<morpho::component_tree<P, AssociativeMap>, V>
    {
    private:
      typedef morpho::component_tree<P, AssociativeMap> tree_t;

    public:
      typedef typename tree_t::vertex_id_t		key_type;
      typedef V						value_type;
      typedef typename std::vector<V>::reference	reference;
      typedef typename std::vector<V>::const_reference  const_reference;

      property_map() = default;

      property_map(const tree_t& tree)
      : m_data(tree._get_data()),
	m_root(tree.get_root_id()),
	m_val(tree.size())
      {
      }

      property_map(const tree_t& tree, const V& init)
      : m_data(tree._get_data()),
	m_root(tree.get_root_id()),
	m_val(tree.size(), init)
      {
      }


      reference operator[] (key_type v)
      {
	mln_precondition(v < m_val.size());
	return m_val[v];
      }

      const_reference operator[] (key_type v) const
      {
	mln_precondition(v < m_val.size());
	return m_val[v];
      }

      reference operator[] (const typename tree_t::node_type& v)
      {
	mln_precondition(v.id() < m_val.size());
	return m_val[v.id()];
      }

      const_reference operator[] (const typename tree_t::node_type& v) const
      {
	mln_precondition(v.id() < m_val.size());
	return m_val[v.id()];
      }



    private:
      const typename tree_t::_data_t*	m_data;
      typename tree_t::vertex_id_t	m_root;
      std::vector<V>			m_val;
    };


}

#endif // ! MLN_MORPHO_COMPONENT_TREE_ATTRIBUTE_MAP_HPP
