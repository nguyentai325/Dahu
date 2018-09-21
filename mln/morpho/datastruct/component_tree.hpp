#ifndef MLN_MORPHO_DATASTRUCT_COMPONENT_TREE_HPP
# define MLN_MORPHO_DATASTRUCT_COMPONENT_TREE_HPP

# include <vector>
# include <memory>
# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/assert.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/iterator/iterator_base.hpp>
# include <mln/core/image/image2d.hpp>

namespace mln
{

  namespace morpho
  {

    /// \brief Data structure that encodes any morphological tree
    /// e.g. mintree, maxtree, tree of shapes. They states the inclusion
    /// of connected components.
    ///
    ///
    /// Implementation details:
    /// A component tree is basically a triplet (N, S, pmap)
    /// where:
    /// + N is a vector of nodes
    /// + S is a vector of points
    /// + pmap is mapping S -> N
    ///
    /// A node is also a triplet (parent, size, first_point_index)
    /// + parent is the parent index (npos for roots)
    /// + size is the number of nodes in the subtree (>= 1)
    /// + first_point_index is the index of the first proper_point in S
    ///
    /// N and S are supposed to be ordered by depth-first search, but it might be too strong
    /// for some application and requires an extra step to produce this ordering, thus
    /// the structure has a tag, that tells if the ordering as been computed. If not, and the
    /// user uses a method that requires this ordering, the ordering process occurs.
    ///
    /// The vector N has an extra sentinel node at the end to ease traversal processes.
    /// This sentinel is at index (npos) and is composed by the triplet
    ///   (parent: npos, prev: root, next: npos, sexts: npos, first_point: S.size())
    /// This sentinal can be considered as the NULL pointer i.e.:
    ///
    /// You can traverse a branch upward:
    /// while (x.id() != NULL) // or x == tree.nend()
    ///   ...
    ///   x = x.parent();

    /// The root node is the triplet (parent:npos, prev: npos, next/nexts: ? size:..., first_point: 0)
    template <class P, class AssociativeMap>
    struct component_tree;

    namespace internal
    {
      // FWD.
      struct component_tree_node;

      template <class P, class AssociativeMap>
      struct component_tree_data;
    }


    namespace internal
    {
      struct component_tree_node
      {
	unsigned    m_parent;	    // The parent node index
	unsigned    m_prev;	    // The previous node index
	unsigned    m_next;         // The next node index (child or next sibling)
	unsigned    m_next_sibling; // The next sibling in the subtree

	//unsigned    m_size;		// The number of nodes in the subtree (>=1)
	unsigned    m_point_index;	// The first point index (in S) of the component
      };



      template <class P, class AssociativeMap>
      struct component_tree_data
      {
	std::vector<component_tree_node>	m_nodes;  // A sorted vector of nodes
	std::vector<P>				m_S; // A sorted vector of points
	AssociativeMap				m_pmap;  // An associative structure P -> node_id_t
	image2d<float_t>   m_Uv;
	image2d<point2d>  m_parent_pixel;
    //image2d<rgb8>  m_U;
	bool					m_pset_ordered = false;
      };
    }



    template <class P, class AssociativeMap>
    struct component_tree
    {
      typedef unsigned		vertex_id_t;
      typedef unsigned		size_type;
      typedef P			point_type;

      component_tree();


      /// \defgroup Site set type and iterator
      /// \{
      struct pset_iterator;
      struct pset_reverse_iterator;
      struct pset_range;
      struct proper_pset_range;

      pset_range pset() const;
      /// \}

      /// \defgroup Node type, iterator & ranges
      /// \{
      struct node_type;
      struct node_iterator;
      struct node_reverse_iterator;
      struct node_range;

      size_type		size() const;
      size_type		realsize() const;
      node_type		get_node(vertex_id_t v) const;
      node_type		get_node_at(const point_type& p) const;
      component_tree	get_subtree(vertex_id_t v) const;
      node_range	nodes(bool ignore_root = false) const;
      node_range	nodes_without_root() const;

      ///  \brief Shrink the nodes vector to free space.
      ///
      /// This is usually used after a filtering operation
      /// to remove unused node from the node vector.
      void              shrink_to_fit();

      /// \}

      /// \defgroup Misc
      /// \{
      struct children_iterator;
      struct children_range;

      static
      constexpr
      vertex_id_t npos()
      {
	return 0; // sentinel index.
      }

      node_type nend() const
      {
	return node_type(this, npos());
      }
      /// \}

      vertex_id_t	get_node_id(const point_type& p) const;
      vertex_id_t	get_root_id() const;
      node_type		get_root() const;

      /// \{
      bool operator== (const component_tree& other) const;
      bool operator!= (const component_tree& other) const;
      /// \}


      /// \defgroup Internal
      /// \{
      typedef morpho::internal::component_tree_data<P, AssociativeMap> _data_t;

      /// \brief Reorder the pset with depth first traversal ordering
      /// Reordering ensures that the first point (canonical element) of the
      /// node stays the same.
      /// FIXME: it actually use the image interface instead of the associative map
      ///        interface. Fix it.
      void		_reorder_pset();

      _data_t*		    _get_data();
      const _data_t*	_get_data() const;
      /// \}

    private:
      component_tree(const _data_t* data, vertex_id_t root);

    private:
      std::shared_ptr<_data_t>  m_data;
      vertex_id_t		m_root;  // The root node index;
    };

    /*************************************************/
    /** Pset range/iterators implementation	    **/
    /*************************************************/
    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::pset_iterator
      : mln::iterator_base<pset_iterator, P, const P&>
    {
      pset_iterator(const _data_t* data, size_type begin, size_type end)
	: m_data (data), m_begin(begin), m_end(end)
      {
      }

      void init()
      {
	m_cur = m_begin;
      }

      void next()
      {
	++m_cur;
      }

      bool finished() const
      {
	return m_cur == m_end;
      }

      const P& dereference() const
      {
	return m_data->m_S[m_cur];
      }

    private:
      const _data_t*	m_data;    // Pointer to component tree data
      size_type		m_begin;   // Start position
      size_type		m_end;	   // Past the end position
      size_type		m_cur;	   // Current position
    };

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::pset_reverse_iterator
      : mln::iterator_base<pset_reverse_iterator, P, const P&>
    {
      pset_reverse_iterator(const _data_t* data, size_type begin, size_type end)
	: m_data (data), m_begin(begin), m_end(end)
      {
      }

      void init()
      {
	m_cur = m_end;
      }

      void next()
      {
	--m_cur;
      }

      bool finished() const
      {
	return m_cur == m_begin;
      }

      const P& dereference() const
      {
	return m_data->m_S[m_cur-1];
      }

    private:
      const _data_t*	m_data;    // Pointer to component tree data
      size_type		m_begin;   // Start position
      size_type		m_end;	   // Past-the-end forward position
      size_type		m_cur;	   // Current position
    };



    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::pset_range
    {
      typedef pset_iterator		iterator;
      typedef pset_iterator		const_iterator;
      typedef pset_reverse_iterator	reverse_iterator;
      typedef pset_reverse_iterator	const_reverse_iterator;

      pset_range(const _data_t* data, vertex_id_t id)
	: m_data(data)
      {
	mln_precondition(id < m_data->m_nodes.size());

	auto x = m_data->m_nodes[id];
	m_begin = x.m_point_index;
	m_end = m_data->m_nodes[x.m_next_sibling].m_point_index;
      }

      iterator iter() const
      {
	return iterator(m_data, m_begin, m_end);
      }

      reverse_iterator riter() const
      {
	return reverse_iterator(m_data, m_begin, m_end);
      }

    private:
      const _data_t*	m_data;    // Pointer to component tree data
      size_type		m_begin;   // Start position
      size_type		m_end;	   // Past-the-end position
    };

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::proper_pset_range
    {
      typedef pset_iterator		iterator;
      typedef pset_iterator		const_iterator;
      typedef pset_reverse_iterator	reverse_iterator;
      typedef pset_reverse_iterator	const_reverse_iterator;

      proper_pset_range(const _data_t* data, vertex_id_t id)
	: m_data(data)
      {
	mln_precondition(id < m_data->m_nodes.size());

	auto x = m_data->m_nodes[id];
	m_begin = x.m_point_index;
	m_end = m_data->m_nodes[x.m_next].m_point_index;
      }

      iterator iter() const
      {
	return iterator(m_data, m_begin, m_end);
      }

      reverse_iterator riter() const
      {
	return reverse_iterator(m_data, m_begin, m_end);
      }

    private:
      const _data_t*	m_data;    // Pointer to component tree data
      size_type		m_begin;   // Start position
      size_type		m_end;	   // Past-the-end position
    };

    /*************************************************/
    /** Children range/iterators implementation	    **/
    /*************************************************/

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::children_iterator
      : mln::iterator_base<children_iterator, node_type, node_type>
    {
      children_iterator() = default;

      children_iterator(const component_tree* tree, vertex_id_t x)
        : m_tree(tree), m_data (tree->_get_data())
      {
        m_begin = m_data->m_nodes[x].m_next;
        m_end = m_data->m_nodes[x].m_next_sibling;
      }

      void init()
      {
        m_cur = m_begin;
      }

      void next()
      {
        m_cur = m_data->m_nodes[m_cur].m_next_sibling;
      }

      bool finished() const
      {
        return m_cur == m_end;
      }

      node_type dereference() const
      {
        return node_type(m_tree, m_cur);
      }

    private:
      const component_tree* m_tree;
      const _data_t*	m_data;    // Pointer to component tree data
      vertex_id_t	m_begin;   // Start position
      vertex_id_t	m_end;	   // Past the end position
      vertex_id_t	m_cur;	   // Current position
    };

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::children_range
    {
      typedef children_iterator iterator;
      typedef children_iterator const_iterator;

      children_range(const component_tree* tree, vertex_id_t id)
        : m_tree(tree), m_node_id(id)
      {
      }

      const_iterator iter() const
      {
        return iterator(m_tree, m_node_id);
      }

    private:
      const component_tree*	m_tree;
      vertex_id_t		m_node_id;	// Vertex id
    };



    /****************************************************/
    /**  Node implementation			        */
    /****************************************************/


    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::node_type
    {
      typedef component_tree<P, AssociativeMap> tree_t;

      node_type()
        : m_tree(NULL)
      {
      }

      node_type(const component_tree* tree, vertex_id_t id)
        : m_tree(tree), m_data (tree->m_data.get()), m_node_id (id)
      {
      }

      bool
      operator== (const node_type& other) const
      {
	return m_node_id == other.m_node_id;
      }

      bool
      operator!= (const node_type& other) const
      {
	return m_node_id != other.m_node_id;
      }

      vertex_id_t id() const
      {
	return m_node_id;
      }

      vertex_id_t get_next_node_id() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_next;
      }

      vertex_id_t get_prev_node_id() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_prev;
      }

      vertex_id_t get_next_sibling_id() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_next_sibling;
      }

      vertex_id_t get_parent_id() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_parent;
      }

      vertex_id_t get_point_index() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_point_index;
      }
      vertex_id_t get_first_child_id() const
      {
	mln_precondition(m_tree);
	if (get_next_node_id() == npos())
	  return npos();
	else if (m_data->m_nodes[get_next_node_id()].m_parent == m_node_id)
	  return get_next_node_id();
	else
	  return npos();
      }


      bool has_child() const
      {
	mln_precondition(m_tree);
	return get_next_node_id() != npos() and
	  m_data->m_nodes[get_next_node_id()].m_parent == m_node_id;
      }

      bool has_next_sibling() const
      {
	mln_precondition(m_tree);
	return get_next_sibling_id() != npos();
      }

      bool is_root() const
      {
	mln_precondition(m_tree);
	return m_node_id == m_tree->m_root;
      }

      node_type	  parent() const
      {
	mln_precondition(m_tree);
	return node_type(m_tree, this->get_parent_id());
      }

      node_type	  next_node() const
      {
	mln_precondition(m_tree);
	return node_type(m_tree, this->get_next_node_id());
      }

      node_type	  prev_node() const
      {
	mln_precondition(m_tree);
	return node_type(m_tree, this->get_prev_node_id());
      }

      node_type	  next_sibling() const
      {
	mln_precondition(m_tree);
	return node_type(m_tree, this->get_next_sibling_id());
      }

      node_type	  first_child() const
      {
	mln_precondition(m_tree);
	return node_type(m_tree, this->get_first_child_id());
      }

      proper_pset_range	  proper_pset() const
      {
	// You must reorder the pset to use this fonction
	// Consider using component_tree::_order_pset()
	mln_precondition(m_tree);
	mln_precondition(m_data->m_pset_ordered);

	return proper_pset_range(m_data, m_node_id);
      }

      pset_range	  pset() const
      {
	// You must reorder the pset to use this fonction
	// Consider using component_tree::_order_pset()
	mln_precondition(m_tree);
	mln_precondition(m_data->m_pset_ordered);

	return pset_range(m_data, m_node_id);
      }

      size_type		  get_first_point_id() const
      {
	mln_precondition(m_tree);
	return m_data->m_nodes[m_node_id].m_point_index;
      }

      P			  first_point() const
      {
	mln_precondition(m_tree);
	return  m_data->m_S[m_data->m_nodes[m_node_id].m_point_index];
      }

      children_range	  children() const
      {
	mln_precondition(m_tree);
	return children_range(m_tree, m_node_id);
      }

      const component_tree&	 tree() const
      {
	mln_precondition(m_tree);
	return *m_tree;
      }

    private:
      const component_tree*	m_tree;    // Pointer to the component tree
      const _data_t*		m_data;    // Pointer to component tree data
      vertex_id_t		m_node_id; // Current vertex position
    };

    /****************************************************/
    /**  Node iterator / range				*/
    /****************************************************/

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::node_iterator
      : iterator_base< node_iterator, node_type, node_type >
    {
      node_iterator() = default;

      node_iterator(const component_tree* tree, vertex_id_t id, bool ignore_root = false)
        : m_tree(tree), m_data(m_tree->m_data.get())
      {
        m_begin = ignore_root ? m_data->m_nodes[id].m_next : id;
        m_end = m_data->m_nodes[id].m_next_sibling;
      }

      void init()
      {
        m_current = m_begin;
      }

      void next()
      {
        m_current = m_data->m_nodes[m_current].m_next;
      }

      bool finished() const
      {
        return m_current == m_end;
      }

      node_type dereference() const
      {
        return node_type(m_tree, m_current);
      }

    private:
      const component_tree* m_tree;
      const _data_t*	    m_data;
      vertex_id_t	m_begin;
      vertex_id_t	m_end;
      vertex_id_t	m_current;
    };

    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::node_reverse_iterator
      : iterator_base< node_reverse_iterator, node_type, node_type >
    {
      node_reverse_iterator() = default;

      node_reverse_iterator(const component_tree* tree, vertex_id_t id, bool ignore_root = false)
        : m_tree(tree), m_data(m_tree->m_data.get())
      {
        m_begin = m_data->m_nodes[ m_data->m_nodes[id].m_next_sibling ].m_prev;
        m_end = ignore_root ? id : m_data->m_nodes[id].m_prev;
      }

      void init()
      {
        m_current = m_begin;
      }

      void next()
      {
        m_current = m_data->m_nodes[m_current].m_prev;
      }

      bool finished() const
      {
        return m_current == m_end;
      }

      node_type dereference() const
      {
        return node_type(m_tree, m_current);
      }

    private:
      const component_tree* m_tree;
      const _data_t*	    m_data;
      vertex_id_t	m_begin;
      vertex_id_t	m_end;
      vertex_id_t	m_current;
    };


    template <class P, class AssociativeMap>
    struct component_tree<P, AssociativeMap>::node_range
    {
      typedef node_iterator		const_iterator;
      typedef node_iterator		iterator;
      typedef node_reverse_iterator	reverse_iterator;
      typedef node_reverse_iterator	const_reverse_iterator;

      node_range(const component_tree* tree, vertex_id_t id, bool ignore_root = false)
        : m_tree(tree), m_node_id(id), m_ignore_root (ignore_root)
      {
      }

      iterator iter() const
      {
        return iterator(m_tree, m_node_id, m_ignore_root);
      }

      reverse_iterator riter() const
      {
        return reverse_iterator(m_tree, m_node_id, m_ignore_root);
      }

      size_type size() const
      {
        return m_tree->size();
      }

      bool has(const node_type& p) const
      {
        if (p.id() == npos())
          return false;

        return (p.tree().get_root_id() == m_node_id); // the same subtree
        // FIXME: if they are not referring the same tree, we should check
        // the current root is in the parenthood of p.
      }

    private:
      const component_tree*	m_tree;
      vertex_id_t		m_node_id;
      bool			m_ignore_root;
    };

    /***************************************/
    /***  Tree method implementation	 ***/
    /***************************************/

    template <class P, class AssociativeMap>
    inline
    component_tree<P, AssociativeMap>::component_tree()
      : m_data(new _data_t()),
        m_root(0)
    {
    }

    template <class P, class AssociativeMap>
    inline
    component_tree<P, AssociativeMap>::component_tree(const _data_t* data, vertex_id_t root)
      : m_data(const_cast<_data_t*>(data)),
	m_root(root)
    {
    }


    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::size_type
    component_tree<P, AssociativeMap>::size() const
    {
      return m_data->m_nodes.size();
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::size_type
    component_tree<P, AssociativeMap>::realsize() const
    {
      size_type n = 1;
      unsigned next = m_data->m_nodes[m_root].m_next;
      while (next != npos()) {
        next = m_data->m_nodes[next].m_next;
        ++n;
      }

      return n;
    }



    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::node_range
    component_tree<P, AssociativeMap>::nodes(bool ignore_root) const
    {
      return node_range(this, m_root, ignore_root);
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::node_range
    component_tree<P, AssociativeMap>::nodes_without_root() const
    {
      return node_range(this, m_root, true);
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::pset_range
    component_tree<P, AssociativeMap>::pset() const
    {
      return pset_range(m_data.get(), m_root);
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::node_type
    component_tree<P, AssociativeMap>::get_node(vertex_id_t v) const
    {
      mln_precondition(v < size());

      return node_type(this, v);
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::vertex_id_t
    component_tree<P, AssociativeMap>::get_node_id(const point_type& p) const
    {
      return m_data->m_pmap[p];
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::node_type
    component_tree<P, AssociativeMap>::get_node_at(const point_type& p) const
    {
      return node_type(this, m_data->m_pmap[p]);
    }

    template <class P, class AssociativeMap>
    inline
    component_tree<P, AssociativeMap>
    component_tree<P, AssociativeMap>::get_subtree(vertex_id_t v) const
    {
      mln_precondition(v < size());

      component_tree x = *this;
      x.m_root = v;
      return x;
    }

    template <class P, class AssociativeMap>
    inline
    const typename component_tree<P, AssociativeMap>::_data_t*
    component_tree<P, AssociativeMap>::_get_data() const
    {
      return m_data.get();
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::_data_t*
    component_tree<P, AssociativeMap>::_get_data()
    {
      return m_data.get();
    }


    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::vertex_id_t
    component_tree<P, AssociativeMap>::get_root_id() const
    {
      return m_root;
    }

    template <class P, class AssociativeMap>
    inline
    typename component_tree<P, AssociativeMap>::node_type
    component_tree<P, AssociativeMap>::get_root() const
    {
      return get_node(m_root);
    }

    template <class P, class AssociativeMap>
    inline
    bool
    component_tree<P, AssociativeMap>::operator== (const component_tree& other) const
    {
      return other.m_data == this->m_data;
    }

    template <class P, class AssociativeMap>
    inline
    bool
    component_tree<P, AssociativeMap>::operator!= (const component_tree& other) const
    {
      return other.m_data != this->m_data;
    }

    template <class P, class AssociativeMap>
    void
    component_tree<P, AssociativeMap>::_reorder_pset()
    {
      mln_entering("mln::morpho::component_tree::_reorder_pset");

      P nullp = value_traits<P>::max();

      std::vector<P> nmap(m_data->m_nodes.size(), nullp);

      mln_ch_value(AssociativeMap, P) next;
      resize(next, m_data->m_pmap);

      mln_reverse_foreach(auto p, m_data->m_S)
	{
	  unsigned k = m_data->m_pmap[p];
	  next[p] = nmap[k];
	  nmap[k] = p;
	}

      unsigned spos = 0;
      mln_foreach(auto x, this->nodes())
	{
	  m_data->m_nodes[x.id()].m_point_index = spos;
	  for (P p = nmap[x.id()]; p != nullp; p = next[p])
	    m_data->m_S[spos++] = p;
	}
      mln_assertion(spos == m_data->m_S.size());

      m_data->m_pset_ordered = true;
      mln_exiting();
    }



    template <class P, class AssociativeMap>
    void
    component_tree<P, AssociativeMap>::shrink_to_fit()
    {
      mln_entering("mln::morpho::component_tree::shrink_to_fit");

      std::vector<vertex_id_t> newidx;
      newidx.resize(m_data->m_nodes.size());

      unsigned n = 1;
      mln_foreach(auto x, this->nodes())
        newidx[x.id()] = n++;
      newidx[0] = 0;

      // Create a new vector of nodes
      // Note: maybe it can be done inplace i.e. in-place permutation ?
      std::vector<internal::component_tree_node> nvec(n);
      {
        int i = 1;
        mln_foreach(auto x, this->nodes())
          {
            nvec[i].m_parent = newidx[x.get_parent_id()];
            nvec[i].m_prev = i-1;
            nvec[i].m_next = i+1;
            nvec[i].m_next_sibling = newidx[x.get_next_sibling_id()];
            nvec[i].m_point_index = x.get_first_point_id();
            i++;
          }

        // Fix last node
        nvec[n-1].m_next = 0;

        // Sentinel
        nvec[0] = {
          0,   // parent -> itself
          n-1, // prev -> lastnode
          0,   // next -> itself
          0,   // next_sibling -> itself
          (unsigned) m_data->m_S.size() // point_index: past the end
        };
      }

      // Update pmap to the new node indexes
      mln_foreach(auto px, m_data->m_pmap.pixels())
        px.val() = newidx[px.val()];

      // Swap buffer
      m_data->m_nodes = std::move(nvec);

      m_root = newidx[m_root];

      mln_exiting();
    }

  }

}

#endif // ! MLN_MORPHO_DATASTRUCT_COMPONENT_TREE_HPP
