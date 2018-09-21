#ifndef MLN_MORPHO_CANVAS_DILATION_LIKE_SPE_HPP
# define MLN_MORPHO_CANVAS_DILATION_LIKE_SPE_HPP

# include <mln/core/win2d.hpp>
# include <mln/core/algorithm/transpose.hpp>

namespace mln
{

  namespace morpho
  {

    namespace canvas
    {

      namespace overload
      {

        // Special case when the SE is a rectangle (seperable)
        // Question ? Should we transpose the image to improve
        // data locality given than inplace transposition
        // might be costly due to cycle detection.
        template <class I, class Compare, class J, class OpTraits>
        typename
        std::enable_if<std::is_same<typename I::domain_type, box2d>::value>::type
        dilation_like(const Image<I>& ima_,
                      const rect2d& nbh,
                      Compare cmp,
                      Image<J>& output,
                      OpTraits __op__)
        {
          box2d r = nbh.dpoints;
          const I& ima = exact(ima_);

          if (r.shape()[0] == 1) {
            impl::dilate_like_0(ima,
                                nbh,
                                cmp,
                                exact(output),
                                __op__,
                                typename image_has_extension<I>::type ());
          } else if (r.shape()[1] == 1) {
            rect2d h0 { box2d {{0, r.pmin[0]}, {1, r.pmax[0]}} };
            mln_concrete(I) tmp = transpose(ima);
            mln_concrete(I) out = imconcretize(tmp);
            morpho::canvas::dilation_like(tmp, h0, cmp, out, __op__);
            transpose(out, output);
          } else {

            rect2d h0 { box2d {{0, r.pmin[0]}, {1, r.pmax[0]}} };
            rect2d h1 { box2d {{0, r.pmin[1]}, {1, r.pmax[1]}} };

            image2d<mln_value(I)> f;
            {
              mln_concrete(I) tmp = imconcretize(ima);
              morpho::canvas::dilation_like(ima, h1, cmp, tmp, __op__);
              f = transpose(tmp);
            }

            {
              mln_concrete(I) tmp = imconcretize(f);
              morpho::canvas::dilation_like(f, h0, cmp, tmp, __op__);
              transpose(tmp, output);
            }
          }
        }

      } // end of namespace mln::morpho::canvas::overload
    } // end of namespace mln::morpho::canvas
  } // end of namespace mln::morpho
} // end of namespace mln

#endif //!MLN_MORPHO_CANVAS_DILATION_LIKE_SPE_HPP
