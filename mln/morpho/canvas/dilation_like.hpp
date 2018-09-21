#ifndef MLN_MORPHO_CANVAS_DILATION_LIKE_HPP
# define MLN_MORPHO_CANVAS_DILATION_LIKE_HPP

# include <mln/core/image/image.hpp>
# include <mln/morpho/se/se.hpp>
# include <mln/kernelv2/kernel.hpp>
# include <mln/core/extension/extension.hpp>
# include <mln/core/extension/fill.hpp>

namespace mln
{

  namespace morpho
  {

    namespace canvas
    {

      //
      // struct dilation_like_operations_traits
      // {
      //   typedef ... support_incremental;
      //   typedef ... aggregate_type;
      //   typedef ... incremental_aggregate_type;
      //   ... zero() constexpr;
      // };

      template <class I, class SE, class Compare, class J, class OpTraits>
      void
      dilation_like(const Image<I>& ima,
                    const StructuringElement<SE>& nbh,
                    Compare cmp,
                    Image<J>& output,
                    OpTraits __op__);

      /******************************************/
      /****          Implementation          ****/
      /******************************************/

      namespace impl
      {

        /// \brief Incremental implementation that enables the dilation to be in
        /// a complexity that does not depend on the SE size.
        /// This specialization is used when:
        /// * The SE is incremental
        /// * The feature has the `untake` method
        template <class I, class SE, class Compare, class J, class OpTraits>
        void
        dilate_like_1(const I& ima, const SE& nbh, Compare cmp, J& out,
                      OpTraits __op__,
                      std::true_type __is_incremental__)
        {
          namespace ker = mln::kernel;
          (void) __is_incremental__;
          (void) __op__;

          ker::Aggregate<typename OpTraits::incremental_aggregate_type> A;
          ker::Point p;
          ker::Neighbor n;
          auto f = ker::make_image_expr<0>(ima);
          auto g = ker::make_image_expr<1>(out);
          auto expr = kernel::declare(g(p) = A (f(n)));
          ker::execute_incremental(expr, nbh);
        }

        /// \brief Basic implementation
        /// This specialization is used when:
        /// * Either the SE is not incremental or nor the feature has the `untake` method.
        template <class I, class SE, class Compare, class J, class OpTraits>
        void
        dilate_like_1(const I& ima, const SE& nbh, Compare cmp, J& out,
                      OpTraits __op__,
                      std::false_type __is_incremental__)
        {
          namespace ker = mln::kernel;
          (void) __is_incremental__;
          (void) __op__;

          ker::Aggregate<typename OpTraits::aggregate_type> A(cmp);
          ker::Point p;
          ker::Neighbor n;
          auto f = ker::make_image_expr<0>(ima);
          auto g = ker::make_image_expr<1>(out);

          auto expr = kernel::declare(g(p) = A (f(n)));
          ker::execute(expr, nbh);
        }

        template <class I, class SE, class Compare, class J, class OpTraits>
        void
        dilate_like_0(const I& ima, const SE& nbh, Compare cmp, J& out,
                      OpTraits __op__,
                      std::true_type __has_extension__)
        {
          (void) __has_extension__;

          mln_value(I) v = __op__.zero();


          using is_se_incremental = typename neighborhood_traits<SE>::is_incremental;
          using is_accu_incremental = typename OpTraits::support_incremental;
          using is_image_incremental = typename iterator_traits<mln_pxter(I)>::has_NL;
          using is_incremental = std::integral_constant<
            bool, is_se_incremental::value and is_accu_incremental::value and
            is_image_incremental::value>;

          if (not is_se_incremental::value)
            mln::trace::warn("Slow because SE is not incremental");
          else if (not is_accu_incremental::value)
            mln::trace::warn("Slow because the accu is not incremental with this image");
          else if (not is_image_incremental::value)
            mln::trace::warn("Slow because the image has no New Line Support");

          if (extension::need_adjust(ima, nbh)) {
            mln::trace::warn("Slow version because input image extension is not wide enough.");
            dilate_like_1(extension::add_value_extension(ima, v), nbh, cmp,
                          out, __op__, is_incremental ());
          } else {
            extension::fill(ima, v);
            dilate_like_1(ima, nbh, cmp, out, __op__, is_incremental ());
          }
        }

        template <class I, class SE, class Compare, class J, class OpTraits>
        void
        dilate_like_0(const I& ima, const SE& nbh, Compare cmp, J& out,
                      OpTraits __op__,
                      std::false_type __has_extension__)
        {
          (void) __has_extension__;
          mln::trace::warn("Slow version because input image has no extension.");


          using is_se_incremental = typename neighborhood_traits<SE>::is_incremental;
          using is_accu_incremental = typename OpTraits::support_incremental;
          using is_image_incremental = typename iterator_traits<mln_pxter(I)>::has_NL;
          using is_incremental = std::integral_constant<
            bool, is_se_incremental::value and is_accu_incremental::value and
            is_image_incremental::value>;

          if (not is_se_incremental::value)
            mln::trace::warn("Slow because SE is not incremental");
          else if (not is_accu_incremental::value)
            mln::trace::warn("Slow because the accu is not incremental with this image");
          else if (not is_image_incremental::value)
            mln::trace::warn("Slow because the image has no New Line Support");

          mln_value(I) v = __op__.zero();
          dilate_like_1(extension::add_value_extension(ima, v), nbh, cmp, out,
                        __op__,
                        is_incremental () );
        }

      } // end of namespace mln::morpho::canvas::impl


      namespace overload
      {
        // Generic implementation
        template <class I, class SE, class Compare, class J, class OpTraits>
        void
        dilation_like(const Image<I>& ima,
                      const StructuringElement<SE>& nbh,
                      Compare cmp,
                      Image<J>& output,
                      OpTraits __op__)
        {
          impl::dilate_like_0(exact(ima),
                              exact(nbh),
                              cmp,
                              exact(output),
                              __op__,
                              typename image_has_extension<I>::type ());
        }
      }

    } // end of namespace mln::morpho::canvas

  } // end of namespace mln::morpho

} // end of namespace mln

# include <mln/morpho/canvas/dilation_like.spe.hpp>

namespace mln
{

  namespace morpho
  {

    namespace canvas
    {

      template <class I, class SE, class Compare, class J, class OpTraits>
      void
      dilation_like(const Image<I>& ima,
                    const StructuringElement<SE>& nbh,
                    Compare cmp,
                    Image<J>& output,
                    OpTraits __op__)
      {
        overload::dilation_like(exact(ima),
                                exact(nbh),
                                cmp,
                                exact(output),
                              __op__);
      }



    } // end of namespace mln::morpho::canvas

  } // end of namespace mln::morpho

} // end of namespace mln



#endif //!MLN_MORPHO_CANVAS_DILATION_LIKE_HPP
