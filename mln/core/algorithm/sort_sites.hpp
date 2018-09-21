#ifndef MLN_CORE_ALGORITHM_SORT_SITES_HPP
# define MLN_CORE_ALGORITHM_SORT_SITES_HPP

# include <mln/core/value/value_traits.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/value/indexer.hpp>
# include <vector>
# include <algorithm>

// FIXME: Speed up HQ version

namespace mln
{


  /// \brief Sort the points of an image.
  ///
  /// \p sort_sites sorts the points of an image. An optional
  /// comparison function can be provided through the parameter \p
  /// cmp.  When two sites \p p and \p q are equivalent (neither
  /// `input(p) < input(q)` nor `input(q) < input(p)` they are ranked
  /// according to the natural point ordering relation (e.g. the scan
  /// order for \p point2d)
  ///
  template <typename I, typename BinaryFunction = std::less<mln_value(I)> >
  std::vector<typename I::site_type>
  sort_sites(const Image<I>& input,  BinaryFunction cmp = BinaryFunction ());


  /****************/
  /* Implem       */
  /****************/

  namespace impl
  {
    template <typename I, typename Compare>
    std::vector<typename I::site_type>
    sort_sites(const I& input, Compare, std::true_type _use_indexer_)
    {
      typedef indexer<mln_value(I), Compare> Indexer;
      typedef typename Indexer::index_type index_t;

      mln_entering("mln::sort_sites (couting sort)");
      (void) _use_indexer_;
      Indexer f;

      typedef mln_value(I) V;
      static constexpr std::size_t nvalues = 1 << value_traits<index_t>::quant;
      unsigned h[nvalues] = {0,};
      {
        mln_viter(v, input);
        mln_forall(v)
          ++h[f(*v)];

        unsigned count = 0;
        index_t i = value_traits<index_t>::min();
        do {
          unsigned tmp = h[i];
          h[i] = count;
          count += tmp;
        } while (i++ < value_traits<index_t>::max());
        assert(count == input.domain().size());
      }

      std::vector<typename I::site_type> v;
      v.resize(input.domain().size());
      {
        mln_pixter(px, input);
        mln_forall(px)
          v[ h[f(px->val())]++ ] = px->point();
      }
      mln_exiting();
      return v;
    }

    template <typename I, typename StrictWeakOrdering>
    std::vector<typename I::site_type>
    sort_sites(const I& input, StrictWeakOrdering cmp, std::false_type _use_indexer_)
    {
      mln_entering("mln::sort_sites (std sort)");

      (void) _use_indexer_;

      std::vector<mln_point(I)> v;
      v.reserve(input.domain().size());
      mln_piter(p, input);
      mln_forall(p)
        v.push_back(*p);

      std::stable_sort(v.begin(), v.end(), [&input, cmp](const mln_point(I)& x, const mln_point(I)& y) {
          return cmp(input(x), input(y)); });

      mln_exiting();
      return v;
    }

  } // end of namespace mln::impl


  template <typename I, typename Compare>
  std::vector<typename I::site_type>
  sort_sites(const Image<I>& input, Compare cmp)
  {
    constexpr bool _is_low_quant = value_traits<mln_value(I)>::quant <= 16;
    constexpr bool _has_indexer  = has_indexer<mln_value(I), Compare>::value;

    std::integral_constant<bool, _is_low_quant and _has_indexer> _use_indexer;
    return impl::sort_sites(exact(input), cmp, _use_indexer);
  }

} // end of namespace mln


#endif // !MLN_CORE_ALGORITHM_SORT_SITES_HPP
