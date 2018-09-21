#ifndef MLN_CORE_CONCEPT_ITERATOR_HPP
# define MLN_CORE_CONCEPT_ITERATOR_HPP

# include <mln/core/concept/object.hpp>
# include <boost/concept_check.hpp>

namespace mln
{

  // Fwd declaration
  template <typename I>
  struct Iterator_;

  template <typename I>
  struct Iterator : Object<I>
  {
    BOOST_CONCEPT_USAGE(Iterator)
    {
      BOOST_CONCEPT_ASSERT((Iterator_<I>));
    }
  };


  // Real concept here.
  template <typename I>
  struct Iterator_
  {
    typedef typename I::value_type    value_type;
    typedef typename I::reference     reference;
    typedef typename I::has_NL        has_NL;

    BOOST_CONCEPT_USAGE(Iterator_)
    {
      // Remove concept because lambda not assignable
      //BOOST_CONCEPT_ASSERT((boost::Assignable<I>));
      //BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<I>));

      reference (I::*method) () const = &I::operator*;
      (void) method;
      void (I::*method1) () = &I::next;
      (void) method1;
      void (I::*method2) () = &I::init;
      (void) method2;
      bool (I::*method3) () const = &I::finished;
      (void) method3;
      I (I::*method4) () const = &I::iter;
      (void) method4;
    }

  };

} // end of namespace mln

#endif //!MLN_CORE_CONCEPT_ITERATOR_HPP
