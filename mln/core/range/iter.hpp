#ifndef MLN_CORE_RANGE_ITER_HPP
# define MLN_CORE_RANGE_ITER_HPP

# include <mln/core/range/range_traits.hpp>
# include <mln/core/iterator/stditerator.hpp>

namespace mln
{

  namespace rng
  {

    template <typename R>
    typename std::conditional<
      is_a<R, Iterator>::value,
      R&,
      typename range_iterator<R>::type>::type
    iter(R& range);


    template <typename R>
    typename range_const_iterator<R>::type
    iter(const R& range);

    template <typename R>
    typename std::conditional<
      is_a<R, Iterator>::value,
      R&,
      typename range_reverse_iterator<R>::type>::type
    riter(R& range);


    template <typename R>
    typename range_const_reverse_iterator<R>::type
    riter(const R& range);


    /*********************/
    /** Implementation  **/
    /*********************/

    namespace impl
    {
      // This is a real MLN range
      template <typename R>
      typename range_iterator<R>::type
      iter(R& range, typename std::enable_if<
	     is_mln_range<R>::value and
	     !is_a<R, Iterator>::value
	     >::type* = NULL)
      {
	return range.iter();
      }

      /// This is a pseudo-range (MLN Iterator)
      template <typename R>
      typename range_iterator<R>::type&
      iter(R& range, typename std::enable_if<
	     is_a<R, Iterator>::value>::type* = NULL)
      {
	return range.iter(); // i.e itself
      }

      // This is a STL range
      template <typename R>
      typename range_iterator<R>::type
      iter(R& range, typename std::enable_if<
	     !is_mln_range<R>::value>::type* = NULL)
      {
	typename range_iterator<R>::type x(range.begin(), range.end());
	return x;
      }

      // This is a real MLN range
      template <typename R>
      typename range_reverse_iterator<R>::type
      riter(R& range, typename std::enable_if<
	     is_mln_range<R>::value and
	     !is_a<R, Iterator>::value
	     >::type* = NULL)
      {
	return range.riter();
      }

      /// This is a pseudo-range (MLN Iterator)
      template <typename R>
      typename range_reverse_iterator<R>::type&
      riter(R& range, typename std::enable_if<
	     is_a<R, Iterator>::value>::type* = NULL)
      {
	return range.riter(); // i.e itself
      }

      // This is a STL range
      template <typename R>
      typename range_reverse_iterator<R>::type
      riter(R& range, typename std::enable_if<
	     !is_mln_range<R>::value>::type* = NULL)
      {
	typename range_reverse_iterator<R>::type x(range.rbegin(), range.rend());
	return x;
      }

    }


    template <typename R>
    typename std::conditional<
      is_a<R, Iterator>::value,
      R&,
      typename range_iterator<R>::type>::type
    iter(R& range)
    {
      return impl::iter(range);
    }

    template <typename R>
    typename range_const_iterator<R>::type
    iter(const R& range)
    {
      return impl::iter(range);
    }

    template <typename R>
    typename std::conditional<
      is_a<R, Iterator>::value,
      R&,
      typename range_reverse_iterator<R>::type>::type
    riter(R& range)
    {
      return impl::riter(range);
    }

    template <typename R>
    typename range_const_reverse_iterator<R>::type
    riter(const R& range)
    {
      return impl::riter(range);
    }


  }

}

#endif // ! MLN_CORE_RANGE_ITER_HPP
