#ifndef MLN_CORE_CONCEPT_CHECK_HPP
# define MLN_CORE_CONCEPT_CHECK_HPP

# include <type_traits>
# include <string>
# include <boost/concept_check.hpp>

# define MLN_CONCEPT_ASSERT_IF(TEST, CONCEPT)    \
  mln::internal::mln_check_if_<TEST, CONCEPT> ();


namespace mln
{

  namespace internal
  {
    template <bool cond, typename T>
    struct mln_check_if_
    {
      mln_check_if_() {}
    };

    template <typename T>
    struct mln_check_if_<true, T>
    {
      mln_check_if_() { BOOST_CONCEPT_ASSERT((T)); }
    };
  }

  template <bool condition>
  using check_t = std::integral_constant<bool, condition>;

  inline
  void check(std::true_type)
  {
  }

  inline
  void check(std::true_type, std::string)
  {
  }

  inline
  void check_false(std::false_type)
  {
  }

  inline
  void check_false(std::false_type, std::string)
  {
  }

  template <typename T>
  T&&
  make_object()
  {
    return static_cast<T&&>(*(typename std::remove_reference<T>::type*)(0));
  }


} // end of namespace mln

#endif //!MLN_CORE_CONCEPT_CHECK_HPP
