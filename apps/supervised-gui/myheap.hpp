#ifndef MYHEAP_HPP
# define MYHEAP_HPP

# include <mln/morpho/component_tree/component_tree.hpp>
# include <vector>

struct myheap
{
private:
  typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;
  typedef mln::property_map<tree_t, float> distance_map_t;
  typedef mln::property_map<tree_t, int> position_map_t;
public:
  myheap(const tree_t& tree,
         const distance_map_t& distancemap,
         position_map_t&       positionmap)
    : m_tree(tree),
      m_dmap(distancemap),
      m_pos(positionmap)
  {
  }


  void push(const tree_t::node_type& p)
  {
    mln_precondition(m_pos[p] == -1);
    m_array.push_back(p.id());
    m_pos[p] = m_array.size() - 1;
    update(p);
  }

  tree_t::node_type pop()
  {
    mln_precondition(not m_array.empty());

    std::swap(m_array.front(), m_array.back());
    tree_t::vertex_id_t p = m_array.back();
    m_array.pop_back();
    m_pos[m_array.front()] = 0;
    m_pos[p] = -1;
    if (not m_array.empty())
      heapify(0);

    return m_tree.get_node(p);
  }

  // H
  void heapify(unsigned i)
  {
    tree_t::vertex_id_t p = m_array[i];
    unsigned sz = m_array.size();
    assert(m_pos[p] != -1);
    while (true)
      {
	unsigned minpos = i;
	if (LCHILD(i) < sz and m_dmap[m_array[minpos]] > m_dmap[m_array[LCHILD(i)]])
	  minpos = LCHILD(i);
	if (RCHILD(i) < sz and m_dmap[m_array[minpos]] > m_dmap[m_array[RCHILD(i)]])
	  minpos = RCHILD(i);
	if (minpos == i)
	  break;
	std::swap(m_array[i], m_array[minpos]);
        m_pos[m_array[i]] = i;
        m_pos[m_array[minpos]] = minpos;
        i = minpos;
      }
  }

  void update(const tree_t::node_type& p)
  {
    int i = m_pos[p];
    float myweight = m_dmap[p];
    mln_assertion(0 <= i and i < (int)m_array.size());
    while (i > 0 and myweight < m_dmap[m_array[PAR(i)]]) {
      m_array[i] = m_array[PAR(i)];
      m_pos[m_array[i]] = i;
      i = PAR(i);
    }
    m_array[i] = p.id();
    m_pos[p] = i;
  }

  bool empty() const { return m_array.empty(); }

  std::size_t size() const { return m_array.size(); }

private:
  static unsigned PAR(unsigned i) { return (i-1) / 2; }
  static unsigned LCHILD(unsigned i) { return 2*i+1; }
  static unsigned RCHILD(unsigned i) { return 2*i+2; }

  const tree_t&                         m_tree;
  const distance_map_t&			m_dmap;
  position_map_t&			m_pos;
  std::vector<tree_t::vertex_id_t>      m_array;
};



#endif // ! MYHEAP_HPP
