#ifndef MLN_CORE_RANGE_CATEGORY_HPP
# define MLN_CORE_RANGE_CATEGORY_HPP

# include <boost/range/category.hpp>

namespace mln {

  /// \brief Tags
  /// \{
  struct single_pass_range_tag {};
  struct forward_range_tag : single_pass_range_tag {};
  struct reversible_range_tag : forward_range_tag {};
  struct bidirectional_range_tag : reversible_range_tag {};
  struct random_access_range_tag : bidirectional_range_tag {};
  /// \}


  template <typename Range>
  struct range_category;

  namespace internal
  {
    template <typename Iterator_tag>
    struct boost_2_mlnrange;

    template <>
    struct boost_2_mlnrange< std::input_iterator_tag >
    {
      typedef single_pass_range_tag type;
    };

    template <>
    struct boost_2_mlnrange< std::forward_iterator_tag >
    {
      typedef forward_range_tag type;
    };

    template <>
    struct boost_2_mlnrange< std::bidirectional_iterator_tag >
    {
      typedef bidirectional_range_tag type;
    };

    template <>
    struct boost_2_mlnrange< std::random_access_iterator_tag >
    {
      typedef random_access_range_tag type;
    };
  }


  namespace internal
  {
    template <typename Range>
    struct has_inner_category_tag
    {
      template <typename R>
      static
      std::true_type test (R*, typename R::category* = NULL) { return std::true_type(); }

      template <typename R>
      static
      std::false_type test (...) { return std::false_type(); }

      typedef decltype( test<Range>((Range*)NULL) ) T;
      static const bool value = T::value;
    };

    template <typename Range, bool use_default>
    struct range_category_helper
    {
      typedef typename internal::boost_2_mlnrange< typename boost::range_category<Range>::type >::type category;
    };

    template <typename Range>
    struct range_category_helper<Range, true>
    {
      typedef typename Range::category category;
    };


  }


  template <typename Range>
  struct range_category
  {
    typedef typename internal::range_category_helper<Range, internal::has_inner_category_tag<Range>::value>::category type;
  };


} // end of namespace mln

#endif //!MLN_CORE_RANGE_CATEGORY_HPP
