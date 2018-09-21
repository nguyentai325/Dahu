#ifndef MAXTREE_ROUTINES_HPP
# define MAXTREE_ROUTINES_HPP


namespace mln
{

  namespace morpho
  {

    namespace internal
    {

      template <typename V>
      point2d
      zfind_repr(const image2d<V>& ima, image2d<point2d>& parent, const point2d& p)
      {
	point2d q = parent(p);
	if (q != p and ima(q) == ima(p))
	  return parent(p) = zfind_repr(ima, parent, q);
	else
	  return p;
      }

      template <typename V>
      unsigned
      zfind_repr(const image2d<V>& ima, image2d<typename image2d<V>::size_type >& parent,
		 typename image2d<V>::size_type p)
      {
	typename image2d<V>::size_type q = parent[p];
	if (q != p and ima[q] == ima[p])
	  return parent[p] = zfind_repr(ima, parent, q);
	else
	  return p;
      }

      // template <typename V>
      // point2d
      // zfind_parent(const image2d<V>& ima, image2d<point2d>& parent, const point2d& p)
      // {
      // 	point2d q = parent(p);
      // 	point2d r = parent(q);
      // 	if (q != r and ima(q) == ima(r))
      // 	  return parent(p) = zfind_repr(ima, parent, r);
      // 	else
      // 	  return q;
      // }

      template <typename V>
      std::size_t
      zfind_parent(const image2d<V>& ima, image2d<typename image2d<V>::size_type>& parent, typename image2d<V>::size_type p)
      {
	return parent[p] = zfind_repr(ima, parent, parent[p]);
      }

      template <typename V>
      std::size_t
      zfind_root(const image2d<V>& ima, image2d<typename image2d<V>::size_type>& parent, typename image2d<V>::size_type p)
      {
	typename image2d<V>::size_type q = parent[p];
	if (q != p)
	  return parent[p] = zfind_root(ima, parent, q);
	else
	  return p;
      }


    }


    template <typename V>
    struct MaxtreeCanonizationAlgorithm
    {
      MaxtreeCanonizationAlgorithm(const image2d<V>& ima,
				   image2d<std::size_t>& parent)
	: m_ima (ima), m_parent (parent)
      {
      }


      void
      operator() (const box2d& domain) const
      {
	image2d<std::size_t> parent = m_parent | domain;
	mln_foreach(auto& p, parent.values())
	  p = internal::zfind_repr(m_ima, m_parent, p);
      }

      const image2d<V>&	m_ima;
      image2d<std::size_t>&   m_parent;

    };

  }

}

#endif // ! MAXTREE_ROUTINES_HPP
