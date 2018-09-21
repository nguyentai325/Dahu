#ifndef MLN_BOX_HH
# define MLN_BOX_HH

# include <mln/core/image/image.hpp>
# include <mln/core/point.hpp>
# include <mln/core/image/internal/nested_loop_iterator.hpp>
# include <mln/core/domain/box_iter.hpp>
# include <tbb/tbb_stddef.h>


namespace mln
{
  // Fwd
  //template <typename T, unsigned dim> struct box_iter;
  //template <typename T, unsigned dim> struct box_rev_iter;


  template <typename T, unsigned dim>
  struct strided_box
  {
    typedef point<T, dim>          point_type;
    typedef point<T, dim>          value_type;

    typedef internal::nested_loop_iterator<
      internal::strided_domain_point_visitor_forward< point<T, dim> >,
      internal::no_op_visitor,
      internal::no_op_visitor,
      internal::point_structure<T, dim>,
      internal::deref_return_point_policy> iterator;

    typedef iterator const_iterator;

    typedef internal::nested_loop_iterator<
      internal::strided_domain_point_visitor_backward< point<T, dim> >,
      internal::no_op_visitor,
      internal::no_op_visitor,
      internal::point_structure<T, dim>,
      internal::deref_return_point_policy> reverse_iterator;

    typedef reverse_iterator const_reverse_iterator;

    strided_box() = default;
    strided_box(const point_type& pmin, const point_type& pmax, const point_type& strides);

    bool	has(const point_type& p) const;
    point_type  shape()	const;
    bool        empty()	const;
    unsigned    size()	const;
    bool	__is_valid() const;

    iterator		iter() const;
    reverse_iterator    riter() const;

    bool		operator== (const strided_box& other) const;
    bool		operator!= (const strided_box& other) const;

    point_type pmin;
    point_type pmax;
    point_type strides;
  };



  //
  // \pre \tparam T must be an integral type
  template <typename T, unsigned dim>
  struct box
  {
    typedef point<T, dim>          point_type;
    typedef point<T, dim>          value_type;

    typedef internal::nested_loop_iterator<
      internal::domain_point_visitor_forward< point<T, dim> >,
      internal::no_op_visitor,
      internal::no_op_visitor,
      internal::point_structure<T, dim>,
      internal::deref_return_point_policy> iterator;

    typedef iterator const_iterator;

    typedef internal::nested_loop_iterator<
      internal::domain_point_visitor_backward< point<T, dim> >,
      internal::no_op_visitor,
      internal::no_op_visitor,
      internal::point_structure<T, dim>,
      internal::deref_return_point_policy> reverse_iterator;


    box() = default;


    constexpr box(const point_type& pmin_, const point_type& pmax_)
      : pmin (pmin_), pmax(pmax_)
    {
    }



    bool has(const point_type& p) const
    {
      mln_precondition(__is_valid());
      for (unsigned i = 0; i < dim; ++i)
	if (p[i] < pmin[i] or p[i] >= pmax[i])
	  return false;
      return true;
    }

    point_type shape() const
    {
      mln_precondition(__is_valid());
      return pmax - pmin;
    }

    bool empty() const
    {
      mln_precondition(__is_valid());
      for (unsigned i = 0; i < dim; ++i)
	if (pmax[i] <= pmin[i])
	  return true;
      return false;
    }

    unsigned size() const
    {
      mln_precondition(__is_valid());
      unsigned sz = 1;
      for (unsigned i = 0; i < dim; ++i)
	sz *= pmax[i] - pmin[i];
      return sz;
    }

    iterator iter() const
    {
      mln_precondition(__is_valid());
      return iterator( internal::point_structure<T, dim> (),
		       internal::make_point_visitor_forward(pmin, pmax),
		       internal::no_op_visitor (),
		       internal::no_op_visitor ());
    }

    reverse_iterator riter() const
    {
      mln_precondition(__is_valid());
      return reverse_iterator( internal::point_structure<T, dim> (),
			       internal::make_point_visitor_backward(pmin, pmax),
			       internal::no_op_visitor (),
			       internal::no_op_visitor ());
    }


#ifndef MLN_DISABLE_PARALLELIZATION

    box(box& r, tbb::split)
    {
      mln_precondition(r.is_divisible());
      pmin = r.pmin;
      pmax = r.pmax;
      r.pmax[0] = r.pmin[0] + (pmax[0] - pmin[0]) / 2 ;
      pmin[0] = r.pmax[0];
    }

#endif

    bool is_divisible() const
    {
      mln_precondition(__is_valid());
      return (pmax[0] - pmin[0]) > 1;
    }

    void join(const box& other)
    {
      mln_precondition(__is_valid());
      mln_precondition(pmax[0] == other.pmin[0]);
      pmax[0] = other.pmax[0];
    }

    bool __is_valid() const
    {
      for (unsigned i = 0; i < dim; ++i)
	if (pmin[i] > pmax[i])
	  return false;
	else if (pmin[i] == pmax[i]) // empty <=> pmin = pmax
	  return pmin == pmax;
      return true;
    }


    bool operator== (const box& other) const
    {
      return pmin == other.pmin and pmax == other.pmax;
    }

    bool operator!= (const box& other) const
    {
      return pmin != other.pmin or pmax != other.pmax;
    }


    point_type pmin;
    point_type pmax;
  };




  typedef box<short, 1> box1d;
  typedef box<float, 1> box1df;
  typedef box<short, 2> box2d;
  typedef box<float, 2> box2df;
  typedef box<short, 3> box3d;

  typedef strided_box<short, 1> sbox1d;
  typedef strided_box<short, 2> sbox2d;
  typedef strided_box<short, 3> sbox3d;



  template <typename T, unsigned dim>
  inline
  std::ostream&
  operator<< (std::ostream& s, const box<T, dim>& b)
  {
    return s << "[" << b.pmin << " ... " << b.pmax << "]";
  }

  template <typename T, unsigned dim>
  inline
  std::ostream&
  operator<< (std::ostream& s, const strided_box<T, dim>& b)
  {
    return s << "[" << b.pmin << " ... " << b.pmax << "..." << b.strides << "]";
  }


  /// Traits

  /// Forward
  template <typename> struct image2d;

  template <>
  struct image_from_domain<box2d>
  {
    template <typename T>
    struct apply
    {
      typedef image2d<T> type;
    };
  };

  template <>
  struct image_from_domain<sbox2d>
  {
    template <typename T>
    struct apply
    {
      typedef image2d<T> type;
    };
  };



  template <typename T, unsigned dim>
  struct grain_box : box<T, dim>
  {
  private:
    typedef box<T, dim> base;

  public:
    grain_box() = default;
    grain_box(const grain_box&) = default;

    grain_box(const box<T, dim>& b, unsigned grain = 1)
      : base(b), m_grain (grain)
    {
      mln_assertion(grain > 0);
    }

    bool is_divisible() const
    {
      mln_precondition(this->__is_valid());
      return (this->pmax[0] - this->pmin[0]) > (int)m_grain;
    }

#ifndef MLN_DISABLE_PARALLELIZATION

    grain_box(grain_box& r, tbb::split)
      : base(r, tbb::split ()), m_grain (r.m_grain)
    {
    }

#endif


  private:
    unsigned m_grain;
  };

  typedef grain_box<short, 1> grain_box1d;
  typedef grain_box<float, 1> grain_box1df;
  typedef grain_box<short, 2> grain_box2d;
  typedef grain_box<float, 2> grain_box2df;
  typedef grain_box<short, 3> grain_box3d;



  /**************************/
  /** Implementation        */
  /**************************/


  template <typename T, unsigned dim>
  strided_box<T, dim>::strided_box(const point_type& pmin_, const point_type& pmax_, const point_type& strides_)
    : pmin (pmin_), strides (strides_)
  {
    mln_precondition(pmin_ <= pmax_);
    for (unsigned i = 0; i < dim; ++i) {
      T q = (pmax_[i] - pmin_[i]) % strides_[i];
      if (q > 0)
	pmax[i] = pmax_[i] - q + strides_[i];
      else
	pmax[i] = pmax_[i];
    }
  }

  template <typename T, unsigned dim>
  inline
  bool
  strided_box<T, dim>::has(const point_type& p) const
  {
    mln_precondition(__is_valid());
    for (unsigned i = 0; i < dim; ++i)
      if (p[i] < pmin[i] or p[i] >= pmax[i] or (p[i] - pmin[i]) % strides[i] != 0)
	return false;
    return true;
  }

  template <typename T, unsigned dim>
  inline
  bool
  strided_box<T, dim>::__is_valid() const
    {
      for (unsigned i = 0; i < dim; ++i)
	if (pmin[i] > pmax[i])
	  return false;
	else if (pmin[i] == pmax[i]) // empty <=> pmin = pmax
	  return pmin == pmax;
      return true;
    }

  template <typename T, unsigned dim>
  inline
  typename strided_box<T, dim>::point_type
  strided_box<T, dim>::shape() const
  {
    point_type shp((pmax - pmin).as_vec() / strides.as_vec());
    return shp;
  }

  template <typename T, unsigned dim>
  inline
  bool
  strided_box<T, dim>::empty() const
  {
    for (unsigned i = 0; i < dim; ++i)
      if (pmax[i] < pmin[i] + strides)
	return true;
    return false;
  }

  template <typename T, unsigned dim>
  inline
  typename strided_box<T, dim>::iterator
  strided_box<T, dim>::iter() const
  {
    return  iterator( internal::point_structure<T, dim> (),
		      internal::make_strided_point_visitor_forward(pmin, pmax, strides),
		      internal::no_op_visitor (),
		      internal::no_op_visitor ());
  }

  template <typename T, unsigned dim>
  inline
  typename strided_box<T, dim>::reverse_iterator
  strided_box<T, dim>::riter() const
  {
    return  reverse_iterator( internal::point_structure<T, dim> (),
			      internal::make_strided_point_visitor_backward(pmin, pmax, strides),
			      internal::no_op_visitor (),
			      internal::no_op_visitor ());

  }

  template <typename T, unsigned dim>
  inline
  bool
  strided_box<T, dim>::operator== (const strided_box& other) const
  {
    return pmin == other.pmin and pmax == other.pmax and strides == other.strides;
  }

  template <typename T, unsigned dim>
  inline
  bool
  strided_box<T, dim>::operator!= (const strided_box& other) const
  {
    return pmin != other.pmin or pmax != other.pmax or strides != other.strides;
  }

}


#endif
