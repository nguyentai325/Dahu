#ifndef MLN_MORPHO_STRUCTURAL_ERODE_HPP
# define MLN_MORPHO_STRUCTURAL_ERODE_HPP

# include <mln/morpho/canvas/dilation_like.hpp>
# include <mln/accu/accumulators/h_infsup.hpp>
# include <mln/accu/accumulators/infsup.hpp>

/// \file

namespace mln
{
  namespace morpho
  {
    namespace structural
    {

      /// \brief Dilation by a structuring element
      /// \param ima Input image
      /// \param se Structuring element
      /// \param cmp Comparaison function
      ///
      template <class I, class SE,
                class Compare = productorder_less<mln_value(I)> >
      mln_concrete(I)
      erode(const Image<I>& ima,
            const StructuringElement<SE>& nbh,
            Compare cmp = Compare ());

      template <class I, class SE, class OutputImage,
                class Compare = productorder_less<mln_value(I)> >
      OutputImage&
      erode(const Image<I>& ima,
            const StructuringElement<SE>& nbh,
            Image<OutputImage>& output,
            Compare cmp = Compare ());

      /******************************************/
      /****          Implementation          ****/
      /******************************************/


      namespace internal
      {

        template <class I, class SE, class Compare>
        struct erode_traits
        {
          using support_incremental = std::integral_constant<
            bool,
            std::is_integral<mln_value(I)>::value and
            std::is_same<Compare, productorder_less<mln_value(I)> >::value and
            (value_traits<mln_value(I)>::quant <= 16)
            >;

          using aggregate_type      = accu::accumulators::inf<mln_value(I), Compare>;
          using incremental_aggregate_type = accu::accumulators::h_inf<mln_value(I)>;


          static
          constexpr
          mln_value(I) zero()
          {
            return value_traits<mln_value(I), Compare>::sup();
          }
        };

      }

      template <class I, class SE, class OutputImage,
                class Compare>
      OutputImage&
      erode(const Image<I>& ima,
            const StructuringElement<SE>& nbh,
            Image<OutputImage>& output,
            Compare cmp)
      {
        mln_entering("mln::morpho::structural::erode");

        morpho::canvas::dilation_like(ima, nbh, cmp, output, internal::erode_traits<I, SE, Compare> ());

        mln_exiting();
        return exact(output);
      }

      template <class I, class SE, class Compare>
      mln_concrete(I)
      erode(const Image<I>& ima,
             const StructuringElement<SE>& nbh,
             Compare cmp)
      {
        mln_concrete(I) out = imconcretize(exact(ima));
        erode(ima, nbh, out, cmp);
        return out;
      }



    } // end of namespace mln::morpho::structural
  } // end of namespace mln::morpho
} // end of namespace mln

#endif //!MLN_MORPHO_STRUCTURAL_ERODE_HPP
