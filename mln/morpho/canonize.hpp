#ifndef CANONIZE_HPP
# define CANONIZE_HPP

# include <mln/morpho/maxtree_routines.hpp>

namespace mln
{
  namespace morpho
  {

    namespace internal
    {
      template <typename V>
      struct Canonizer
      {
	typedef typename image2d<V>::size_type size_type;

	Canonizer(const image2d<V>& ima, image2d<size_type>& parent, size_type* S)
	  : m_ima(ima), m_parent(parent), m_out (S)
	{
	  resize(m_dejavu, m_ima).init(false);
	  //m_back = S + ima.domain().size();
	}


	void canonize(size_type p)
	{
	  m_dejavu[p] = true;
	  size_type q = m_parent[p];
	  if (not m_dejavu[q])
	    canonize(q);
	  if (m_ima[q] == m_ima[m_parent[q]])
	    m_parent[p] = m_parent[q];
	  *(m_out++) = p;
	}


	// void canonize(size_type p)
	// {
	//   size_type q = internal::zfind_repr(m_ima, m_parent, p);
	//   size_type q0 = q;

	//   int i = 0;
	//   while (!m_dejavu[q]) {
	//     m_dejavu[q] = true;
	//     m_parent[q] = internal::zfind_repr(m_ima, m_parent, m_parent[q]);
	//     q = m_parent[q];
	//     ++i;
	//   }

	//   if (p != q0) {
	//     m_parent[p] = q0;
	//     *(--m_back) = p;
	//   }

	//   for (int j = i-1; j >= 0; --j) {
	//     m_front[j] = q0;
	//     q0 = m_parent[q0];
	//   }
	//   m_front += i;
	//   m_dejavu[p] = true;
	// }


	const image2d<V>&	m_ima;
	image2d<size_type>&	m_parent;
	size_type*		m_out;
	// size_type*		m_front;
	// size_type*		m_back;
	image2d<bool>		m_dejavu;
      };

    }

    template <typename V>
    void
    canonize(const image2d<V>& ima,
	     image2d< typename image2d<V>::size_type >& parent,
	     typename image2d<V>::size_type*	   S)
    {
      internal::Canonizer<V> o(ima, parent, S);

      mln_pixter(p, ima);
      mln_forall(p)
	if (not o.m_dejavu[p->index()])
	  o.canonize(p->index());
    }


    template <typename V>
    void
    canonize(const image2d<V>& ima,
	     const std::vector<typename image2d<V>::size_type>& S,
	     image2d<typename image2d<V>::size_type>& parent)
    {
      typedef typename image2d<V>::size_type size_type;
      for(size_type p: S)
	{
	  size_type q = parent[p];
	  if (ima[q] == ima[parent[q]])
	    parent[p] = parent[q];
	}
    }

  }
}

#endif // ! CANONIZE_HPP
