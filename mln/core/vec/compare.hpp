#ifndef MLN_CORE_VEC_COMPARE_HPP
# define MLN_CORE_VEC_COMPARE_HPP

# include <mln/core/value/value_traits.hpp>

namespace mln
{

  /// \brief Object function for a comparison between
  /// two vectorial types using the lexicographical order.
  ///
  /// \f[
  /// a < b \leftrightarrow \exists i s.t. a_i < b_i
  /// \text{ and } \forall j < i \, \neg b_j < a_j
  /// \f]
  template <typename U, typename V>
  struct lexicographicalorder_less;

  template <typename U, typename V>
  bool
  veclex_isless(const U&, const V&);

  template <typename U, typename V>
  bool
  veclex_islessequal(const U&, const V&);

  template <typename U, typename V>
  bool
  veclex_isgreater(const U&, const V&);

  template <typename U, typename V>
  bool
  veclex_isgreaterequal(const U&, const V&);


  /// \brief Object function for a comparison between
  /// two vectorial types using the product order.
  ///
  /// \f[
  /// a \leq b \leftrightarrow \forall i \; a_i \leq b_i
  /// \f]
  /// Strict relation are deduced by reflexive reduction i.e
  ///
  /// \f[
  /// a < b \leftrightarrow a \leq b and a \ne b
  /// \f]
  template <typename U, typename V>
  struct productorder_less;

  template <typename U, typename V>
  struct productorder_less_equal;

  template <typename U, typename V = U>
  struct productorder_greater;

  template <typename U, typename V = U>
  struct productorder_greater_equal;

#if 0 // Facade only

  template <typename U, typename V>
  bool
  vecprod_isless(const U&, const V&);

  template <typename U, typename V>
  bool
  vecprod_islessequal(const U&, const V&);

  template <typename U, typename V>
  bool
  vecprod_isgreater(const U&, const V&);

  template <typename U, typename V>
  bool
  vecprod_isgreaterequal(const U&, const V&);

#endif

  /*********************/
  /**   Traits        **/
  /*********************/

  template <typename U, unsigned dim, typename tag>
  struct value_traits< internal::vec_base<U, dim, tag>,
                       std::less<internal::vec_base<U, dim, tag> > >
  {
  private:
    typedef internal::vec_base<U, dim, tag> Vec;

    static
    void __check()
    {
      static_assert(internal::vec_base_traits<tag>::is_less_than_comparable,
                    "This type does not support comparison by default."
                    "You must explicity provide the comparison function.");
    }

  public:
    static constexpr unsigned quant = value_traits<U>::quant;
    static constexpr unsigned ndim = dim;

    static constexpr Vec min()
    {
      return __check(), Vec( value_traits<U>::min() );
    }

    static constexpr Vec max()
    {
      return __check(), Vec( value_traits<U>::max() );
    }

    static constexpr Vec inf()
    {
      return __check(), min();
    }

    static constexpr Vec sup()
    {
      return __check(), max();
    }

  };


  template <typename U, unsigned dim, typename tag>
  struct value_traits< internal::vec_base<U, dim, tag>,
		       lexicographicalorder_less< internal::vec_base<U, dim, tag> >
		       >
  {
  private:
    typedef internal::vec_base<U, dim, tag> Vec;

  public:
    static constexpr unsigned quant = value_traits<U>::quant;
    static constexpr unsigned ndim = dim;

    static constexpr Vec min()
    {
      return Vec( value_traits<U>::min() );
    }

    static constexpr Vec max()
    {
      return Vec( value_traits<U>::max() );
    }

    static constexpr Vec inf()
    {
      return min();
    }

    static constexpr Vec sup()
    {
      return max();
    }
  };


  template <typename U, unsigned dim, typename tag>
  struct value_traits< internal::vec_base<U, dim, tag>,
		       productorder_less< internal::vec_base<U, dim, tag> >
		       >
  {
  private:
    typedef internal::vec_base<U, dim, tag> Vec;

  public:
    static constexpr unsigned quant = value_traits<U>::quant;
    static constexpr unsigned ndim = dim;

    static constexpr Vec min()
    {
      return Vec( value_traits<U>::min() );
    }

    static constexpr Vec max()
    {
      return Vec( value_traits<U>::max() );
    }

    static constexpr Vec inf()
    {
      return min();
    }

    static constexpr Vec sup()
    {
      return max();
    }
  };


  /*************************/
  /**  inf/sup overloads  **/
  /*************************/

  template <typename U, unsigned dim, typename tag>
  internal::vec_base<U, dim, tag>
  inf(const internal::vec_base<U, dim, tag>& u,
      const internal::vec_base<U, dim, tag>& v,
      productorder_less< internal::vec_base<U, dim, tag> >)
  {
    internal::vec_base<U, dim, tag> res;
    for (unsigned i = 0; i < dim; ++i)
      res[i] = std::min(u[i], v[i]);
    return res;
  }

  template <typename U, unsigned dim, typename tag>
  internal::vec_base<U, dim, tag>
  inf(const internal::vec_base<U, dim, tag>& u,
      const internal::vec_base<U, dim, tag>& v)
  {
	return inf(u, v, productorder_less< internal::vec_base<U, dim, tag> > ());
  }


  template <typename U, unsigned dim, typename tag>
  internal::vec_base<U, dim, tag>
  sup(const internal::vec_base<U, dim, tag>& u,
      const internal::vec_base<U, dim, tag>& v,
      productorder_less< internal::vec_base<U, dim, tag> >)
  {
    internal::vec_base<U, dim, tag> res;
    for (unsigned i = 0; i < dim; ++i)
      res[i] = std::max(u[i], v[i]);
    return res;
  }

  template <typename U, unsigned dim, typename tag>
  internal::vec_base<U, dim, tag>
  sup(const internal::vec_base<U, dim, tag>& u,
      const internal::vec_base<U, dim, tag>& v)
  {
	return sup(u, v, productorder_less< internal::vec_base<U, dim, tag> > ());
  }

  /*********************/
  /** Implementation  **/
  /*********************/

  template <typename U, typename V, unsigned dim, typename tag>
  bool
  veclex_isless(const internal::vec_base<U, dim, tag>& u,
		const internal::vec_base<V, dim, tag>& v)
  {
    for (unsigned i = 0; i < dim; ++i)
      {
	if (u[i] < v[i])
	  return true;
	else if (v[i] < u[i])
	  return false;
      }
    return false;
  }

  template <typename U, typename V, unsigned dim, typename tag>
  bool
  veclex_islessequal(const internal::vec_base<U, dim, tag>& u,
		     const internal::vec_base<V, dim, tag>& v)
  {
    for (unsigned i = 0; i < dim; ++i)
      {
	if (u[i] < v[i])
	  return true;
	else if (v[i] < u[i])
	  return false;
      }
    return true;
  }

  template <typename U, typename V, unsigned dim, typename tag>
  bool
  veclex_isgreater(const internal::vec_base<U, dim, tag>& u,
		   const internal::vec_base<V, dim, tag>& v)
  {
    return veclex_isless(v, u);
  }


  template <typename U, typename V, unsigned dim, typename tag>
  bool
  veclex_isgreaterequal(const internal::vec_base<U, dim, tag>& u,
			const internal::vec_base<V, dim, tag>& v)
  {
    return veclex_islessequal(v, u);
  }

  template <typename U, typename V, unsigned dim, typename tag>
  bool
  vecprod_isless(const internal::vec_base<U, dim, tag>& u,
		 const internal::vec_base<V, dim, tag>& v)
  {
    bool res = false;
    for (unsigned i = 0; i < dim; ++i) {
      if (v[i] < u[i]) // neg: u[i] <= v[i]
	return false;
      res |= (u[i] != v[i]); // reflexive reduction
    }
    return res;
  }

  template <typename U, typename V>
  typename std::enable_if< std::is_arithmetic<U>::value and
                           std::is_arithmetic<V>::value, bool>::type
  vecprod_isless(const U& u, const V& v)
  {
    return u < v;
  }



  template <typename U, typename V, unsigned dim, typename tag>
  bool
  vecprod_islessequal(const internal::vec_base<U, dim, tag>& u,
		      const internal::vec_base<V, dim, tag>& v)
  {
    for (unsigned i = 0; i < dim; ++i)
      if (v[i] < u[i]) // beg: u[i] <= v[i]
	return false;
    return true;
  }

  template <typename U, typename V>
  typename std::enable_if< std::is_arithmetic<U>::value and
                           std::is_arithmetic<V>::value, bool>::type
  vecprod_islessequal(const U& u, const V& v)
  {
    return u <= v;
  }


  template <typename U, typename V, unsigned dim, typename tag>
  bool
  vecprod_isgreater(const internal::vec_base<U, dim, tag>& u,
		    const internal::vec_base<V, dim, tag>& v)
  {
    return vecprod_isless(v, u);
  }


  template <typename U, typename V, unsigned dim, typename tag>
  bool
  vecprod_isgreaterequal(const internal::vec_base<U, dim, tag>& u,
			 const internal::vec_base<V, dim, tag>& v)
  {
    return vecprod_islessequal(v, u);
  }


  template <typename U, typename V>
  struct lexicographicalorder_less
  {
    bool
    operator() (const U& u, const V& v) const
    {
      return u < v;
    }
  };


  template <typename U, typename V,
	    unsigned dim, typename tag>
  struct lexicographicalorder_less< internal::vec_base<U, dim, tag>, internal::vec_base<V, dim, tag> >
  {
    bool
    operator() (const internal::vec_base<U, dim, tag>& u,
		const internal::vec_base<V, dim, tag>& v) const
    {
      return veclex_isless(u,v);
    }
  };

  template <typename U, typename V>
  struct productorder_less
  {
    bool
    operator() (const U& u, const V& v) const
    {
      return vecprod_isless(u, v);
    }
  };

  template <typename U, typename V>
  struct productorder_less_equal
  {
    bool
    operator() (const U& u, const U& v) const
    {
      return vecprod_islessequal(u, v);
    }
  };

  template <typename U, typename V,
	    unsigned dim, typename tag>
  struct productorder_greater< internal::vec_base<U, dim, tag>, internal::vec_base<V, dim, tag> >
  {
    bool
    operator() (const internal::vec_base<U, dim, tag>& u,
		const internal::vec_base<V, dim, tag>& v) const
    {
      return vecprod_isgreater(u, v);
    }
  };

  template <typename U, typename V,
	    unsigned dim, typename tag>
  struct productorder_greater_equal< internal::vec_base<U, dim, tag>, internal::vec_base<V, dim, tag> >
  {
    bool
    operator() (const internal::vec_base<U, dim, tag>& u,
		const internal::vec_base<V, dim, tag>& v) const
    {
      return vecprod_isgreaterequal(u, v);
    }
  };

}

#endif // ! MLN_CORE_VEC_COMPARE_HPP
