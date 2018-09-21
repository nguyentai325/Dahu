#ifndef MLN_MORPHO_TOS_CTOS_HPP
# define MLN_MORPHO_TOS_CTOS_HPP

# include <mln/core/image/image.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/tos/impl/ctos_serial.hpp>
# include <mln/morpho/tos/impl/ctos_parallel.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      template <typename Compare>
      struct equiv;
    }


    /// \brief Compute the tree of shapes
    ///
    /// Compute the tree of shapes my immerging the image in Khalimsky grid.
    /// The user is responsible to provide an image with a constant border, since the
    /// algorithm starts flooding at point ima.domain().pmin. The user may also want
    /// to perform a median interpolation to ensure correct properties on TOS.
    ///
    /// \return
    ///
    template <typename I,
              typename Neighborhood,
              typename Compare = std::less<mln_value(I)>,
              typename Equiv = internal::equiv<Compare>,
              bool use_priority = false>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima,
         const Neighborhood& nbh,
         mln_point(I) pmin,
         const Compare& cmp,
         const Equiv& eq);


    template <typename I, typename Neighborhood, typename Compare = std::less<mln_value(I)>,
              typename Equiv = internal::equiv<Compare> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp, const Equiv& eq);

    template <typename I, typename Neighborhood, typename Compare = std::less<mln_value(I)> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp = Compare () );

    template <typename I, typename Neighborhood, typename Compare = std::less<mln_value(I)> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_pinf(const Image<I>& ima, const Neighborhood& nbh, mln_point(I) pmin,
              const Compare& cmp = Compare ());

    /***********************************************/
    /* Same as before but using priority proagation */
    /***********************************************/

    template <typename I,
              typename Neighborhood,
              typename Compare = std::less<mln_value(I)>,
              typename Equiv = internal::equiv<Compare> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, mln_point(I) pmin, const Compare& cmp, const Equiv& eq);

    template <typename I, typename Neighborhood, typename Compare = std::less<mln_value(I)>, typename Equiv = internal::equiv<Compare> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp, const Equiv& eq);

    template <typename I, typename Neighborhood, typename Compare = std::less<mln_value(I)> >
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp = Compare () );




    /********************/
    /** Implementation **/
    /********************/


    namespace internal
    {

      template <typename Compare>
      struct equiv
      {
        equiv(const Compare& cmp) :
        m_cmp (cmp)
        {
        }

        template <typename T>
        bool operator () (const T& x, const T& y) const
        {
          return !m_cmp(x,y) and !m_cmp(y,x);
        }

      private:
        Compare m_cmp;
      };

    }

    template <typename I,
              typename Neighborhood,
              typename Compare,
              typename Equiv,
              bool use_priority>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima,
         const Neighborhood& nbh,
         mln_point(I) pmin,
         const Compare& cmp,
         const Equiv& eq)
    {
      return impl::parallel::cToS(ima, nbh, pmin, cmp, eq);
    }


    template <typename I,
              typename Neighborhood,
              typename Compare,
              typename Equiv>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp, const Equiv& equiv)
    {
      mln_point(I) pmin = exact(ima).domain().pmin;
      return impl::parallel::cToS(ima, nbh, pmin, cmp, equiv);
    }



    template <typename I, typename Neighborhood, typename Compare>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp)
    {
      mln_point(I) pmin = exact(ima).domain().pmin;
      //std::cout << "pmin   "  << pmin  << std::endl;
      return impl::parallel::cToS(ima, nbh, pmin, cmp, internal::equiv<Compare> (cmp));
    }

    template <typename I, typename Neighborhood, typename Compare>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_pinf(const Image<I>& ima, const Neighborhood& nbh, mln_point(I) pmin,
              const Compare& cmp)
    {
      return impl::parallel::cToS(ima, nbh, pmin, cmp, internal::equiv<Compare> (cmp));
    }


    template <typename I,
              typename Neighborhood,
              typename Compare,
              typename Equiv>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, mln_point(I) pmin, const Compare& cmp, const Equiv& equiv)
    {
      return impl::parallel::cToS<I, Neighborhood, Compare, Equiv, true>(ima, nbh, pmin, cmp, equiv);
    }

    template <typename I,
              typename Neighborhood,
              typename Compare,
              typename Equiv>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp, const Equiv& equiv)
    {
      mln_point(I) pmin = exact(ima).domain().pmin;
      return impl::parallel::cToS<I, Neighborhood, Compare, Equiv, true>(ima, nbh, pmin, cmp, equiv);
    }

    template <typename I, typename Neighborhood, typename Compare>
    morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
    cToS_priority(const Image<I>& ima, const Neighborhood& nbh, const Compare& cmp)
    {
      mln_point(I) pmin = exact(ima).domain().pmin;
      return impl::parallel::cToS<I, Neighborhood, Compare, internal::equiv<Compare>, true>(ima, nbh, pmin, cmp, internal::equiv<Compare> (cmp));
    }

  }

}

#endif // ! CTOS_HPP
