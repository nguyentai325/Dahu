#ifndef MLN_MORPHO_FILTERING_HPP
# define MLN_MORPHO_FILTERING_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/colors.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      struct Yes
      {
	template <typename T>
	constexpr bool operator() (T) const
	{
	  return true;
	}
      };
    }


    /// \brief Compute area attribute on a component tree (i.e ToS, min/max tree).
    ///
    /// \param ima Image on which the tree has been computed
    /// \param parent Parent image
    /// \param S S vector
    /// \param is_face_2 Filtering function
    /// \param Compare Strict Weak Ordering, relation used for tree computation
    ///
    /// \return 
    ///
    template <typename V, class is_face_2_t = internal::Yes, typename Compare = std::less<V> >
    image2d<unsigned>
    area_compute(const image2d<V>& ima,
		 const image2d<typename image2d<V>::size_type>& parent,
		 const std::vector<typename image2d<V>::size_type>& S,
		 const is_face_2_t& is_face2 = is_face_2_t (),
		 const Compare& cmp = Compare ());


    // template <typename T, typename V,  class is_face_2_t = internal::Yes, typename Compare = std::less<V>, typename Sum >
    // image2d< T >
    // area_filter(const image2d< T >& ima,
    // 		const image2d<V>& K,
    // 		const image2d<typename image2d<V>::size_type>& parent,
    // 		const std::vector<typename image2d<V>::size_type>& S,
    // 		unsigned fvalue,
    // 		const is_face_2_t& is_face2 = is_face_2_t (),
    // 		const Compare& cmp = Compare ());



    /*************************/
    /** Implementation      **/
    /*************************/

    template <typename V, class is_face_2_t = internal::Yes, typename Compare = std::less<V> >
    image2d<unsigned>
    area_compute(const image2d<V>& ima,
		 const image2d<typename image2d<V>::size_type>& parent,
		 const std::vector<typename image2d<V>::size_type>& S,
		 const is_face_2_t& is_face2 = is_face_2_t (),
		 const Compare& cmp = Compare ())
    {
      typedef typename image2d<V>::size_type size_type;
      image2d<unsigned> area;
      resize(area, ima).init(0);

      point2d q = ima.point_at_index(S[0]);
      area[S[0]] += is_face2(q);
      for (int i = S.size()-1; i > 0; --i) // all except root
	{
	  size_type p = S[i];

	  point2d q = ima.point_at_index(p);
	  if (is_face2(q))
	    area[p] += 1;
	  area[parent[p]] += area[p];
	}

      // propagate to canonical element to non-canonical
      for (size_type p: S)
	{
	  V a = ima[p], b = ima[parent[p]];
	  if (!cmp(a, b) and !cmp(b,a)) // a == b, same node
	    area[p] = area[parent[p]];
	}

      return area;
    }


    template <typename T, typename V, class is_face_2_t = internal::Yes, typename Compare = std::less<V>,
	      typename SumType = rgb<std::size_t> >
    image2d<T>
    area_filter(const image2d<T>& ima,
		const image2d<V>& K,
		const image2d<typename image2d<V>::size_type>& parent,
		const std::vector<typename image2d<V>::size_type>& S,
		unsigned lambda,
		const is_face_2_t& is_face2 = is_face_2_t (),
		const SumType& = SumType(),
		const Compare& cmp = Compare ())
    {
      typedef typename image2d<V>::size_type size_type;

      //mln_precondition(ima.domain() == K.domain());
      //mln_precondition(ima.domain() == parent.domain());

      SumType dval = SumType ();
      image2d<SumType> acc;
      image2d<unsigned> area;
      image2d<unsigned> count;
      image2d<T> out;

      resize(acc, ima).init(dval);
      resize(area, ima).init(0);
      resize(count, ima).init(0);
      resize(out, ima);

      point2d q = ima.point_at_index(S[0]);
      if (is_face2(q)) {
	area[S[0]] = 1;
	acc[S[0]] = ima[S[0]];
      }

      for (int i = S.size()-1; i > 0; --i) // all except root
	{
	  size_type p = S[i];
	  size_type rp = K[p] == K[parent[p]] ? parent[p] : p;
	  if (is_face2(K.point_at_index(p))) { // 2-face
	    area[p] += 1;
	    count[rp] += 1;
	    acc[rp] += ima[p];
	  }
	  area[parent[p]] += area[p];
	  //acc[parent[p]] += acc[p];
	}

      // filter + propagate to canonical element to non-canonical
      std::cout << "Area: " << area[S[0]] << std::endl;
      std::cout << lambda << std::endl;

      // image2d<bool> isleaf;
      // resize(isleaf, ima, ima.border(), false);
      for (size_type p: S)
	{
	  size_type q = parent[p];
	  V a = K[p], b = K[q];
	  if (!cmp(a, b) and !cmp(b,a))
	    {// a == b, same node
	      area[p] = area[q];
	      acc[p] = acc[q];
	    }

	  if (area[p] < lambda) {
	    if (area[q] < lambda)
	      out[p] = out[q];
	    else {
	      while (count[q] == 0)
		q = parent[q];
	      out[p] = acc[q] / count[q];
	      //	      isleaf[q] = true;
	    }
	  } else {
	    out[p] = ima[p];
	  }
	}

      // for (size_type p: S)
      // 	{
      // 	  unsigned rq = K[p] == K[parent[p]] ? parent[p] : p;
      // 	  if (isleaf[rq])
      // 	    out[p] = acc[rq] / count[rq];
      // 	}

      return out;
    }


    template <typename T, typename V, typename is_face_2_t = internal::Yes, typename Compare = std::less<V> >
    image2d< rgb<T> >
    area_filterr(const image2d< rgb<T> >& ima,
		 const image2d<V>& K,
		 const image2d<typename image2d<V>::size_type>& parent,
		 const std::vector<typename image2d<V>::size_type>& S,
		 unsigned lambda,
		 const is_face_2_t& is_face2 = is_face_2_t (),
		 const Compare& cmp = Compare ())
    {
      typedef rgb<std::size_t> SumType;
      typedef typename image2d<V>::size_type size_type;

      //mln_precondition(ima.domain() == K.domain());
      //mln_precondition(ima.domain() == parent.domain());

      SumType dval = SumType ();
      image2d<unsigned> area;
      image2d< rgb<T> > out;


      resize(area, ima, ima.border(), 0);
      resize(out, ima);

      point2d q = ima.point_at_index(S[0]);
      area[S[0]] += (is_face2(q));

      for (int i = S.size()-1; i > 0; --i) // all except root
	{
	  size_type p = S[i];
	  point2d q = K.point_at_index(p);
	  area[p] += is_face2(q);
	  area[parent[p]] += area[p];
	}

      // filter + propagate to canonical element to non-canonical
      std::cout << "Area: " << area[S[0]] << std::endl;
      std::cout << lambda << std::endl;

      out[S[0]] = ima[S[0]];
      for (size_type p: S)
	{
	  size_type q = parent[p];
	  V a = K[p], b = K[q];
	  if ((!cmp(a, b) and !cmp(b,a)))
	    area[p] = area[q];

	  if (area[p] > lambda)
	    out[p] = ima[p];
	  else
	    out[p] = out[q];

	}

      return out;
    }

  }
}


#endif // ! MLN_MORPHO_FILTERING_HPP
