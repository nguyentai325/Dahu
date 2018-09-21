#ifndef MLN_MORPHO_STRUCTURAL_GRADIENT_HPP
# define MLN_MORPHO_STRUCTURAL_GRADIENT_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/algorithm/transform.hpp>
# include <mln/core/math_ops.hpp>
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

      /// \brief Compute the Beucher gradient.
      ///
      /// The beucher gradient is defined as:
      /// \f[
      ///  |\Nabla u(p)| = \left| \bigvee_{q \in \mathcal{N}(p)} f(q) -
      ///    |\bigwedge_{q \in \mathcal{N}(p)} f(q) \right|
      /// \f]
      ///
      /// \param[in] input        Input image.
      /// \param[in] se   Neighborhood/SE/Window to look around.
      /// \param[in] cmp  (optional) Comparaison function. The method internally does an
      ///                 unqualified call to `inf(x, y, cmp)` and `sup(x, y, cmp)`. Default
      ///                 is the product-order.
      /// \param[in] norm (optional) Norm function used in \f$|x - y|\f$. The default is
      ///                  the ℓ₂-norm.
      /// \param[out] out (optional) Output image to write in.
      ///
      /// FIXME: This can be done in a single pass using kernels but you
      /// we have to consider:
      /// * an incrementable accumlator
      /// * an extension that is mirrorizable
      template <class I, class SE,
                class Compare = productorder_less<mln_value(I)>,
                class NormFunction = mln::functional::l2norm_t<mln_value(I)> >
      mln_ch_value(I, typename std::result_of<NormFunction(mln_value(I))>::type)
        gradient(const Image<I>& input,
                 const StructuringElement<SE>& se,
                 Compare cmp = Compare (),
                 NormFunction norm = NormFunction ());

      template <class I, class SE, class Compare, class NormFunction, class O>
      void
      gradient(const Image<I>& input,
               const StructuringElement<SE>& se,
               Compare cmp,
               NormFunction norm,
               Image<O>& out);


      /// \brief Compute the morphological external gradient.
      ///
      /// The beucher external gradient is defined as:
      /// \f[
      ///  |\Nabla u(p)| = \left| \bigvee_{q \in \mathcal{N}(p)} f(q) - f(p) \right|
      /// \f]
      ///
      /// \param[in] input        Input image.
      /// \param[in] se   Neighborhood/SE/Window to look around.
      /// \param[in] cmp  (optional) Comparaison function. The method internally does an
      ///                 unqualified call to `inf(x, y, cmp)` and `sup(x, y, cmp)`. Default
      ///                 is the product-order.
      /// \param[in] norm (optional) Norm function used in \f$|x - y|\f$. The default is
      ///                  the ℓ₂-norm.
      /// \param[out] out (optional) Output image to write in.
      ///
      /// FIXME: This can be done in a single pass using kernels but you
      /// we have to consider:
      /// * an incrementable accumlator
      /// * an extension that is mirrorizable
      template <class I, class SE,
                class Compare = productorder_less<mln_value(I)>,
                class NormFunction = mln::functional::l2norm_t<mln_value(I)> >
      mln_ch_value(I, typename std::result_of<NormFunction(mln_value(I))>::type)
        external_gradient(const Image<I>& input,
                 const StructuringElement<SE>& se,
                 Compare cmp = Compare (),
                 NormFunction norm = NormFunction ());

      template <class I, class SE, class Compare, class NormFunction, class O>
      void
      external_gradient(const Image<I>& input,
                        const StructuringElement<SE>& se,
                        Compare cmp,
                        NormFunction norm,
                        Image<O>& out);

      /*************************/
      /***  Implementation   ***/
      /*************************/

      namespace impl
      {

        // Version non-fast dilate - erode
        template <class I, class SE, class Compare, class Norm, class J>
        void
        gradient(const I& ima, const SE& nbh, Compare cmp, Norm norm, J& out)
        {
          auto d = morpho::structural::dilate(ima, nbh, cmp);
          auto e = morpho::structural::erode(ima, nbh, cmp);

          transform(d-e, norm, out);
        }

        // Version non-fast dilate - erode
        template <class I, class SE, class Compare, class Norm, class J>
        void
        external_gradient(const I& ima, const SE& nbh, Compare cmp, Norm norm, J& out)
        {
          auto d = morpho::structural::dilate(ima, nbh, cmp);

          transform(d-ima, norm, out);
        }


      }

      template <class I, class SE, class Compare, class Norm, class O>
      void
      gradient(const Image<I>& ima_,
               const StructuringElement<SE>& se_,
               Compare cmp,
               Norm norm,
               Image<O>& output)
      {
        mln::trace::entering("mln::morpho::gradient");

        const I& ima = exact(ima_);
        const SE& se = exact(se_);
        O& out = exact(output);
        mln::morpho::structural::impl::gradient(ima, se, cmp, norm, out);

        mln::trace::exiting();
      }


      template <class I, class SE, class Compare, class Norm>
      mln_ch_value(I, typename std::result_of<Norm(mln_value(I))>::type)
        gradient(const Image<I>& ima_,
                 const StructuringElement<SE>& se,
                 Compare cmp, Norm norm)
      {
        const I& ima = exact(ima_);

        typedef mln_value(I) V;
        typedef typename std::result_of<Norm(V)>::type result_type;

        mln_ch_value(I, result_type) out = imchvalue<result_type> (ima);
        mln::morpho::structural::gradient(ima, se, cmp, norm, out);

        return out;
      }

      template <class I, class SE, class Compare, class Norm, class O>
      void
      external_gradient(const Image<I>& ima_,
                        const StructuringElement<SE>& se_,
                        Compare cmp,
                        Norm norm,
                        Image<O>& output)
      {
        mln::trace::entering("mln::morpho::external_gradient");

        const I& ima = exact(ima_);
        const SE& se = exact(se_);
        O& out = exact(output);
        mln::morpho::structural::impl::external_gradient(ima, se, cmp, norm, out);

        mln::trace::exiting();
      }


      template <class I, class SE, class Compare, class Norm>
      mln_ch_value(I, typename std::result_of<Norm(mln_value(I))>::type)
        external_gradient(const Image<I>& ima_,
                          const StructuringElement<SE>& se,
                          Compare cmp, Norm norm)
      {
        const I& ima = exact(ima_);

        typedef mln_value(I) V;
        typedef typename std::result_of<Norm(V)>::type result_type;

        mln_ch_value(I, result_type) out = imchvalue<result_type> (ima);
        mln::morpho::structural::external_gradient(ima, se, cmp, norm, out);

        return out;
      }


    }

  }

}

#endif // !MLN_MORPHO_STRUCTURAL_GRADIENT_HPP
