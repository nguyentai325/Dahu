#ifndef MLN_CORE_OPS_HPP
# define MLN_CORE_OPS_HPP

# include <functional>
# include <algorithm>

/**
* \file
* \brief Defines the fundamental logical, arithmetic and relational operators:
* \see mln::logical_not
* \see mln::logical_and
* \see mln::logical_or
* \see mln::negate
* \see mln::add
* \see mln::substract
* \see mln::multiples
* \see mln::devides
* \see mln::modulo
* \see mln::equals_to
* \see mln::not_equals_to
* \see mln::greater_than
* \see mln::less_than
* \see mln::greater_equal
* \see mln::less_equal
* \see mln::getter
* \see mln::inf
* \see mln::sup
* \see mln::min
* \see mln::max
*/


namespace mln
{
  using std::negate;
  using std::logical_not;

  /**
  * \struct add
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed addition arithmetic operation.
  *
  * \struct substract
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed substraction arithmetic operation.
  *
  * \struct multiplies
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed multiplication arithmetic operation.
  *
  * \struct devides
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed division arithmetic operation.
  *
  * \struct modulo
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed division arithmetic operation.
  *
  * \struct equals_to
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed equality comparison operation.
  *
  * \struct not_equals_to
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed inequality comparison operation.
  *
  * \struct less_than
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed <em>lesser than</em> comparison operation.
  *
  * \struct greater_than
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed <em>greater than</em> comparison operation.
  *
  * \struct less_equal
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed <em>lesser equal</em> comparison operation.
  *
  * \struct greater_equal
  * \tparam U
  * \tparam V
  * \brief A function object for the mixed <em>greater equal</em> comparison operation.
  */

  template <typename U, typename V = U>
  struct add : std::binary_function<U, V, typename std::common_type<U,V>::type>
  {
    typename std::common_type<U,V>::type
    operator() (const U& x, const V& y) const
    {
      return x+y;
    }
  };

  template <typename U, typename V = U>
  struct substract : std::binary_function<U, V, typename std::common_type<U,V>::type>
  {
    typename std::common_type<U,V>::type
    operator() (const U& x, const V& y) const
    {
      return x-y;
    }
  };

  template <typename U, typename V = U>
  struct multiplies : std::binary_function<U, V, typename std::common_type<U,V>::type>
  {
    typename std::common_type<U,V>::type
    operator() (const U& x, const V& y) const
    {
      return x*y;
    }
  };

  template <typename U, typename V = U>
  struct devides : std::binary_function<U, V, typename std::common_type<U,V>::type>
  {
    typename std::common_type<U,V>::type
    operator() (const U& x, const V& y) const
    {
      return x/y;
    }
  };

  template <typename U, typename V = U>
  struct modulo : std::binary_function<U, V, typename std::common_type<U,V>::type>
  {
    typename std::common_type<U,V>::type
    operator() (const U& x, const V& y) const
    {
      return x%y;
    }
  };


  template <typename U, typename V = U>
  struct equal_to : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x == y;
    }
  };

  template <typename U, typename V = U>
  struct not_equal_to : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x != y;
    }
  };


  template <typename U, typename V = U>
  struct less_than : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x < y;
    }
  };


  template <typename U, typename V = U>
  struct greater_than : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x > y;
    }
  };

  template <typename U, typename V = U>
  struct less_equal : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x <= y;
    }
  };

  template <typename U, typename V = U>
  struct greater_equal : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x >= y;
    }
  };

  template <typename U, typename V = U>
  struct logical_and : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x && y;
    }
  };

  template <typename U, typename V = U>
  struct logical_or : std::binary_function<U, V, bool>
  {
    bool
    operator() (const U& x, const V& y) const
    {
      return x || y;
    }
  };

}

// FIXME: get must be imported in the mln namespace to prevent
// a dependant name lookup that would force the lookup at
// definition instead of the instanciation.
namespace mln
{
  using std::get;

  template <size_t N>
  struct getter
  {
    template <class C>
    auto
    operator() (C&& obj) const
      -> decltype(get<N>(std::forward<C>(obj)))
    {
      using namespace std;
      return get<N>(std::forward<C>(obj));
    }
  };

}

  /****************************/
  /**  Relational operations **/
  /****************************/
namespace mln
{
  using std::min;
  using std::max;


  template <typename U>
  const U&
  inf(const U& x, const U& y)
  {
    return std::min(x, y);
  }

  template <typename U>
  const U&
  sup(const U& x, const U& y)
  {
    return std::max(x, y);
  }

  template <typename U, class Compare>
  const U&
  inf(const U& x, const U& y, Compare cmp)
  {
    return std::min(x, y, cmp);
  }

  template <typename U, class Compare>
  const U&
  sup(const U& x, const U& y, Compare cmp)
  {
    return std::max(x, y, cmp);
  }

  namespace functional
  {

    template <class U, class V = U>
    struct max_t : std::binary_function<U, V, typename std::common_type<U,V>::type >
    {
      typedef typename std::common_type<U,V>::type R;

      R
      operator() (const U& x, const V& y) const
      {
        return std::max<R>(x,y);
      }
    };

    template <class U, class V = U>
    struct min_t : std::binary_function<U, V, typename std::common_type<U,V>::type >
    {
      typedef typename std::common_type<U,V>::type R;

      R
      operator() (const U& x, const V& y) const
      {
        return std::min<R>(x,y);
      }
    };

    template <class U, class V = U>
    struct inf_t : std::binary_function<U, V, typename std::common_type<U,V>::type >
    {
      typedef typename std::common_type<U,V>::type R;

      R
      operator() (const U& x, const V& y) const
      {
        return inf<R>(x,y);
      }
    };

    template <class U, class V = U>
    struct sup_t : std::binary_function<U, V, typename std::common_type<U,V>::type >
    {
      typedef typename std::common_type<U,V>::type R;

      R
      operator() (const U& x, const V& y) const
      {
        return sup<R>(x,y);
      }
    };

  }

  /*****************************/
  /** Aggregation operations  **/
  /*****************************/


} // end of namespace mln

#endif //!MLN_CORE_OPS_HPP
