#ifndef MLN_MORPHO_STRUCTURAL_CLOSING_HPP
# define MLN_MORPHO_STRUCTURAL_CLOSING_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/algorithm/transform.hpp>
# include <mln/core/trace.hpp>
# include <mln/morpho/se/se.hpp>
# include <mln/morpho/structural/erode.hpp>
# include <mln/morpho/structural/dilate.hpp>

/// \file

namespace mln
{

  namespace morpho
  {

    namespace structural
    {

      /// \brief Compute the morphological closing.
      ///
      /// \[
      /// \gamma(f) = \varepsilon(\delta(f))
      /// \]
      ///
      /// \param[in] input        Input image.
      /// \param[in] se   Neighborhood/SE/Window to look around.
      /// \param[in] cmp  (optional) Comparaison function. The method internally does an
      ///                 unqualified call to `inf(x, y, cmp)` and `sup(x, y, cmp)`. Default
      ///                 is the product-order.
      /// \param[out] out (optional) Output image to write in.
      ///
      template <class I, class SE,
                class Compare = productorder_less<mln_value(I)> >
      mln_concrete(I)
      closing(const Image<I>& input,
              const StructuringElement<SE>& se,
              Compare cmp = Compare ());

      template <class I, class SE, class Compare, class O>
      O&
      closing(const Image<I>& input,
               const StructuringElement<SE>& se,
               Compare cmp,
               Image<O>& out);

      /*************************/
      /***  Implementation   ***/
      /*************************/

      namespace impl
      {

        // Version non-fast dilate - erode
        template <class I, class SE, class Compare, class J>
        void
        closing(const I& ima, const SE& nbh, Compare cmp, J& out)
        {
          auto d = morpho::structural::dilate(ima, nbh, cmp);
          morpho::structural::erode(d, nbh, out, cmp);
        }

      }

      template <class I, class SE, class Compare, class O>
      O&
      closing(const Image<I>& ima_,
              const StructuringElement<SE>& se_,
              Compare cmp,
              Image<O>& output)
      {
        mln::trace::entering("mln::morpho::closing");

        const I& ima = exact(ima_);
        const SE& se = exact(se_);
        O& out = exact(output);
        mln::morpho::structural::impl::closing(ima, se, cmp, out);

        mln::trace::exiting();
        return out;
      }


      template <class I, class SE, class Compare>
      mln_concrete(I)
      closing(const Image<I>& ima_,
              const StructuringElement<SE>& se,
              Compare cmp)
      {
        const I& ima = exact(ima_);

        mln_concrete(I) out = imconcretize(ima);
        mln::morpho::structural::closing(ima, se, cmp, out);

        return out;
      }

    }

  }

}

#endif // !MLN_MORPHO_STRUCTURAL_CLOSING_HPP
