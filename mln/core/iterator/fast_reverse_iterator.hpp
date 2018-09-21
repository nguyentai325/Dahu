// (C) Copyright David Abrahams 2002.
// (C) Copyright Jeremy Siek    2002.
// (C) Copyright Thomas Witt    2002.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MLN_CORE_ITERATOR_FAST_REVERSE_ITERATOR_HPP
# define MLN_CORE_ITERATOR_FAST_REVERSE_ITERATOR_HPP


#include <boost/iterator.hpp>
#include <boost/utility.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

namespace mln
{

  /// This reverse iterator enables
  /// to write fast reverse iterator that does not create temporaries.
  /// As a counterpart, the base class iterator must have a valid "before begin"
  /// sementic
  template <class Iterator>
  class fast_reverse_iterator
    : public boost::iterator_adaptor< fast_reverse_iterator<Iterator>, Iterator >
  {
    typedef boost::iterator_adaptor< fast_reverse_iterator<Iterator>, Iterator > super_t;

    friend class boost::iterator_core_access;

   public:
    fast_reverse_iterator() {}

    explicit fast_reverse_iterator(Iterator x)
      : super_t(--x) {}

    template<class OtherIterator>
    fast_reverse_iterator(fast_reverse_iterator<OtherIterator> const& r,
                          typename std::enable_if<std::is_convertible<OtherIterator, Iterator>::value>* = 0)
      : super_t(r.base())
    {}

   private:
      typename super_t::reference dereference() const { return *this->base_reference(); }

      void increment() { --this->base_reference(); }
      void decrement() { ++this->base_reference(); }

      void advance(typename super_t::difference_type n)
      {
          this->base_reference() += -n;
      }

      template <class OtherIterator>
      typename super_t::difference_type
      distance_to(fast_reverse_iterator<OtherIterator> const& y) const
      {
        return this->base_reference() - y.base();
      }
  };

} // namespace mln


#endif //! MLN_CORE_ITERATOR_FAST_REVERSE_ITERATOR_HPP
