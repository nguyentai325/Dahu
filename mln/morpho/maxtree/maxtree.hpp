#ifndef MLN_MORPHO_MAXTREE_MAXTREE_HPP
# define MLN_MORPHO_MAXTREE_MAXTREE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/morpho/maxtree/maxtree_queue.hpp>

namespace mln
{

  namespace morpho
  {
    ///
    /// \brief Compute a maxtree as a component tree where nodes
    /// contain indexes instead of points.
    ///
    /// \param ima The image
    /// \param nbh The neighborhood
    /// \param cmp A strict weak order on values
    template <typename I, typename Im, typename Pa, typename N, typename StrictWeakOrdering = std::less<mln_value(I)> >
    component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    maxtree_indexes(const Image<I>& ima, const Image<Im>& ima1, const Image<Pa> & ima2, const Neighborhood<N>& nbh,
                    StrictWeakOrdering cmp = StrictWeakOrdering());

    ///
    /// \brief Compute a mintree as a component tree where nodes
    /// contain indexes instead of points.
    ///
    /// \param ima The image
    /// \param nbh The neighborhood
    /// \param cmp A strict weak order on values
    template <typename I, typename Im, typename Pa, typename N>
    component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    mintree_indexes(const Image<I>& ima, const Image<Im>& ima1, const Image<Pa>& ima2,const Neighborhood<N>& nbh);


    /*****************************/
    /**  Implementation		**/
    /*****************************/

    template <typename I, typename Im, typename Pa, typename N, typename StrictWeakOrdering>
    component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    maxtree_indexes(const Image<I>& ima, const Image<Im>& ima1, const Image<Pa> & ima2, const Neighborhood<N>& nbh, StrictWeakOrdering cmp)
    {
      mln_entering("mln::morpho::maxtree_indexes");
      auto res = impl::maxtree_queue_indexes(exact(ima), exact(ima1), exact(ima2), exact(nbh), cmp);
      mln_exiting();
      return res;
    }

    template <typename I,typename Im, typename Pa ,typename N>
    component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    mintree_indexes(const Image<I>& ima, const Image<Im>& ima1, const Image<Pa>& ima2, const Neighborhood<N>& nbh)
    {
      mln_entering("mln::morpho::mintree_indexes");
      auto res = impl::maxtree_queue_indexes(exact(ima), exact(ima1), exact(ima2), exact(nbh), std::greater<mln_value(I)> ());
      mln_exiting();
      return res;
    }



  }

}

#endif // ! MLN_MORPHO_MAXTREE_MAXTREE_HPP
