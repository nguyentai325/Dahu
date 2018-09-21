#ifndef MLN_CORE_RANGE_ADAPTOR_REVERSED_HPP
# define MLN_CORE_RANGE_ADAPTOR_REVERSED_HPP

# include <mln/core/range/iterator_range.hpp>
# include <mln/core/range/range_reverse_iterator.hpp>
# include <mln/core/range/category.hpp>
# include <mln/core/range/rbegin.hpp>
# include <mln/core/range/rend.hpp>

namespace mln {
  namespace range {
    namespace adaptor {

      template <typename ReversibleRange>
      struct reversed_range
        : public iterator_range<typename range_reverse_iterator<ReversibleRange>::type>
      {
        // FIXME: add concept check
      private:
        typedef iterator_range<typename range_reverse_iterator<ReversibleRange>::type> base;

      public:
        typedef typename range_category<ReversibleRange>::type category;

        explicit reversed_range(ReversibleRange& r)
          : base(rbegin(r), rend(r)) // rbegin, rend found by ADL
        {
        }

      };

      struct reverse_t {};

    } // end of namespace mln::range::adaptor
  } // end of namespace mln::range

    /// \brief operator| overload for reversed adaptor
    template <typename ReversibleRange>
    range::adaptor::reversed_range<ReversibleRange>
    operator | (ReversibleRange& r, range::adaptor::reverse_t)
    {
      return range::adaptor::reversed_range<ReversibleRange>(r);
    }

    inline
    constexpr range::adaptor::reverse_t reversed()
    {
      return range::adaptor::reverse_t();
    }

    template <typename ReversibleRange>
    inline
    range::adaptor::reversed_range<ReversibleRange>
    reversed(ReversibleRange& r)
    {
      return range::adaptor::reversed_range<ReversibleRange>(r);
    }

} // end of namespace mln




#endif //!MLN_CORE_RANGE_ADAPTOR_REVERSED_HPP
