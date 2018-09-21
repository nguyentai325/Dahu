#ifndef RANGE_TRAITS_HPP
# define RANGE_TRAITS_HPP

# include <type_traits>
# include <mln/core/concept/iterator.hpp>
# include <mln/core/iterator/stditerator.hpp>

namespace mln
{
  struct input_range_tag {};
  struct output_range_tag {};
  struct forward_range_tag : input_range_tag {} ;
  struct bidirectional_range_tag : forward_range_tag {};
  struct random_access_range_tag : bidirectional_range_tag {};



  template <typename Range>
  struct range_traits;

  template <typename Range>
  struct is_mln_range;

  template <typename Range, bool mln_range = is_mln_range<Range>::value>
  struct range_iterator;

  template <typename Range, bool mln_range = is_mln_range<Range>::value>
  struct range_const_iterator;

  template <typename Range, bool mln_range = is_mln_range<Range>::value>
  struct range_reverse_iterator;

  template <typename Range, bool mln_range = is_mln_range<Range>::value>
  struct range_const_reverse_iterator;

  template <typename Range>
  struct range_value;

  template <typename Range>
  struct range_reference;

  /**************************/
  /** Implementation       **/
  /**************************/

  namespace impl
  {
    template <typename R, typename dummy = void>
    struct has_iter_member : std::false_type
    {
    };

    template <typename R>
    struct has_iter_member<R, decltype(std::declval<R>().iter(), (void)0)>
      : std::true_type
    {
    };
  }

  template <typename R>
  struct is_mln_range :
    mln::is_a<typename R::iterator, Iterator>
  {
  };

  template <typename R>
  struct range_iterator<R, true>
  {
    typedef typename R::iterator			type;
  };

  template <typename R>
  struct range_iterator<R, false>
  {
    typedef stditerator<typename R::iterator>		type;
  };

  template <typename R>
  struct range_iterator<const R, true>
  {
    typedef typename R::const_iterator			type;
  };

  template <typename R>
  struct range_iterator<const R, false>
  {
    typedef stditerator<typename R::const_iterator>	type;
  };

  template <typename R>
  struct range_const_iterator<R, true>
  {
    typedef typename R::const_iterator			type;
  };

  template <typename R>
  struct range_const_iterator<R, false>
  {
    typedef stditerator<typename R::const_iterator>	type;
  };


  template <typename R>
  struct range_reverse_iterator<R, true>
  {
    typedef typename R::reverse_iterator			type;
  };

  template <typename R>
  struct range_reverse_iterator<R, false>
  {
    typedef stditerator<typename R::reverse_iterator>		type;
  };

  template <typename R>
  struct range_reverse_iterator<const R, true>
  {
    typedef typename R::const_reverse_iterator			type;
  };

  template <typename R>
  struct range_reverse_iterator<const R, false>
  {
    typedef stditerator<typename R::const_reverse_iterator>	type;
  };

  template <typename R>
  struct range_const_reverse_iterator<R, true>
  {
    typedef typename R::const_reverse_iterator			type;
  };

  template <typename R>
  struct range_const_reverse_iterator<R, false>
  {
    typedef stditerator<typename R::const_reverse_iterator>	type;
  };


  template <typename R>
  struct range_value
  {
    typedef typename range_iterator<R>::type::value_type type;
  };

  template <typename R>
  struct range_reference
  {
    typedef typename range_iterator<R>::type::reference	type;
  };

  namespace internal
  {

    template <typename R, bool mln_range>
    struct range_traits;

    template <typename R>
    struct range_traits<R, true>
    {
      typedef typename R::category category;
    };

    template <class iterator_tag>
    struct iterator_to_range_trait;

    template <>
    struct iterator_to_range_trait<std::input_iterator_tag>
    {
      typedef input_range_tag category;
    };

    template <>
    struct iterator_to_range_trait<std::output_iterator_tag>
    {
      typedef output_range_tag category;
    };

    template <>
    struct iterator_to_range_trait<std::forward_iterator_tag>
    {
      typedef forward_range_tag category;
    };

    template <>
    struct iterator_to_range_trait<std::bidirectional_iterator_tag>
    {
      typedef bidirectional_range_tag category;
    };

    template <>
    struct iterator_to_range_trait<std::random_access_iterator_tag>
    {
      typedef random_access_range_tag category;
    };

    template <typename R>
    struct range_traits<R, false>
    {
      typedef typename iterator_to_range_trait<
	typename std::iterator_traits<typename R::iterator>::type
	>::category category;
    };
  }


  template <typename R>
  struct range_traits : internal::range_traits<R, is_mln_range<R>::value>
  {
  };

}

#endif // ! RANGE_TRAITS_HPP
