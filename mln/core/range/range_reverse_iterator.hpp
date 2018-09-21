#ifndef MLN_CORE_RANGE_RANGE_REVERSE_ITERATOR_HPP
# define MLN_CORE_RANGE_RANGE_REVERSE_ITERATOR_HPP

namespace mln {

  template <typename ReversibleRange>
  struct range_reverse_iterator
  {
    typedef typename ReversibleRange::reverse_iterator type;
  };

} // end of namespace mln

#endif //!MLN_CORE_RANGE_RANGE_REVERSE_ITERATOR_HPP
