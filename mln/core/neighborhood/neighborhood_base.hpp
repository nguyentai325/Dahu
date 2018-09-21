#ifndef MLN_CORE_NEIGHBORHOOD_NEIGHBORHOOD_BASE_HPP
# define MLN_CORE_NEIGHBORHOOD_NEIGHBORHOOD_BASE_HPP

# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/core/concept/iterator.hpp>
# include <mln/core/concept/pixel.hpp>
# include <boost/utility/enable_if.hpp>

namespace mln
{
  namespace internal
  {

    template <typename Nbh, typename X>
    struct nbh_bind_point
    {
      typedef decltype(std::declval<Nbh*>()->__bind_point(std::declval<X&>())) type;
    };

    template <typename Nbh, typename X>
    struct nbh_process_point
    {
      typedef decltype(std::declval<Nbh*>()->__process_point(std::declval<X>())) type;
    };

    template <typename Nbh, typename X>
    struct nbh_bind_pixel
    {
      typedef decltype(std::declval<Nbh*>()->__bind_pixel(std::declval<X&>())) type;
    };

    template <typename Nbh, typename X,
	      bool is_pixel = is_a<typename X::value_type, Pixel>::value>
    struct nbh_bind_iterator
    {
      typedef decltype(std::declval<Nbh*>()->__bind_point_iterator(std::declval<X&>())) type;
    };

    template <typename Nbh, typename X>
    struct nbh_bind_iterator<Nbh, X, true>
    {
      typedef decltype(std::declval<Nbh*>()->__bind_pixel_iterator(std::declval<X&>())) type;
    };

  }


  template <class Nbh, class tag>
  struct neighborhood_base : Neighborhood<Nbh>
  {

  public:
    typedef tag                 category;
    typedef std::false_type     is_incremental;

    /// \brief Overload if x is a point lvalue or rvalue
    template <typename X>
    typename boost::lazy_enable_if_c<
      not mln::is_a<X, Iterator>::value and
      not mln::is_a<X, Pixel>::value,
      internal::nbh_bind_point<Nbh, X> >::type
    operator() (X& x) const
    {
      return exact(this)->__bind_point(x);
    }


  /// \brief Overload if x is a point rvalue
  template <typename X>
  typename boost::lazy_enable_if_c<
      not std::is_lvalue_reference<X>::value and
      not mln::is_a<X, Iterator>::value and
      not mln::is_a<X, Pixel>::value,
      internal::nbh_process_point<Nbh, X> >::type
  operator() (X&& x) const
  {
    return exact(this)->__process_point(std::forward<X>(x));
  }

    /// \brief Overload if x is a pixel lvalue
    template <typename X>
    typename boost::lazy_enable_if_c<
      ! mln::is_a<X, Iterator>::value and
        mln::is_a<X, Pixel>::value,
      internal::nbh_bind_pixel<Nbh, X> >::type
    operator() (X& x) const
    {
      return exact(this)->__bind_pixel(x);
    }


    /// \brief Overload if x is a point iterator
    template <typename X>
    typename boost::lazy_enable_if_c<
      !mln::is_a<typename X::value_type, Pixel>::value,
      internal::nbh_bind_iterator<Nbh, X, false>
      >::type
    __dispatch_iterator(X& x) const
    {
      return exact(this)->__bind_point_iterator(x);
    }

    /// \brief Overload if x is a pixel iterator
    template <typename X>
    typename boost::lazy_enable_if_c<
      mln::is_a<typename X::value_type, Pixel>::value,
      internal::nbh_bind_iterator<Nbh, X, true>
      >::type
    __dispatch_iterator (X& x) const
    {
      return exact(this)->__bind_pixel_iterator(x);
    }

    /// \brief Overload if x is an iterator
    template <typename X>
    typename boost::lazy_enable_if_c<
      mln::is_a<X, Iterator>::value,
      internal::nbh_bind_iterator<Nbh, X> >::type
    operator() (X& x) const
    {
      return this->__dispatch_iterator(x);
    }

  };

}

#endif // ! MLN_CORE_NEIGHBORHOOD_NEIGHBORHOOD_BASE_HPP


