#ifndef MLN_CORE_ITERATOR_ITERATOR_TRAITS_HPP
# define MLN_CORE_ITERATOR_ITERATOR_TRAITS_HPP

# include <iterator>

namespace mln
{

  template <typename I>
  struct iterator_traits
  {
    typedef typename I::value_type value_type;
    typedef typename I::reference  reference;

    typedef typename I::has_NL     has_NL;
  };


  template <typename T>
  struct iterator_traits<T*> : std::iterator_traits<T*>
  {
  };

  template <typename T>
  struct iterator_traits<const T*> : std::iterator_traits<const T*>
  {
  };

}
#endif // ! MLN_CORE_ITERATOR_ITERATOR_TRAITS_HPP
