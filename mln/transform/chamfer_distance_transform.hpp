#ifndef MLN_TRANSFORM_CHAMFER_DISTANCE_TRANSFORM_HPP
# define MLN_TRANSFORM_CHAMFER_DISTANCE_TRANSFORM_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace transform
  {

    template <class I, class N, class OutputImage>
    void
    chamfer_distance_transform(const Image<I>& f,
                               const Neighborhood<N>& nbh,
                               OutputImage&& out);


    template <class I, class N, typename DistanceType = int>
    mln_ch_value(I, DistanceType)
    chamfer_distance_transform(const Image<I>& f,
                               const Neighborhood<N>& nbh);


    /**********************************/
    /**   Implementation            ***/
    /**********************************/

    template <class I, class N, class OutputImage>
    void
    chamfer_distance_transform(const Image<I>& f_,
                               const Neighborhood<N>& nbh_,
                               OutputImage&& out)
    {
      static_assert( std::is_convertible<mln_value(I), bool>::value,
                     "Input value type must be convertible to bool");

      mln_entering("mln::transform::chamfer_distance_transform");

      typedef mln_value(OutputImage) distance_t;

      const I& f = exact(f_);
      const N& nbh = exact(nbh_);


      extension::fill(out, 0);

      // Forward scan
      {
        mln_pixter(pxin, pxout, f, out);
        mln_iter(q, nbh(pxout));

        mln_forall(pxin, pxout) {
          if (pxin->val()) {
            distance_t vmin = value_traits<distance_t>::max();
            mln_forall(q)
              if (q->index() < pxout->index())
                vmin = std::min<distance_t>(vmin, q->val() + 1);
            pxout->val() = vmin;
          } else {
            pxout->val() = 0;
          }
        }
      }

      // Backward
      {
        mln_riter(pxin, f.pixels());
        mln_riter(pxout, out.pixels());
        mln_iter(q, nbh(pxout));

        mln_forall(pxin, pxout)
          if (pxin->val())
            {
              distance_t vmin = pxout->val();
              mln_forall(q)
                if (q->index() > pxout->index())
                  vmin = std::min<distance_t>(vmin, q->val()+1);
              pxout->val() = vmin;
            }
      }

      mln_exiting();
    }

    template <class I, class N, typename DistanceType>
    mln_ch_value(I, DistanceType)
    chamfer_distance_transform(const Image<I>& f_,
                               const Neighborhood<N>& nbh_)
    {
      const I& f = exact(f_);
      const N& nbh = exact(nbh_);

      mln_ch_value(I, DistanceType) out = imchvalue<DistanceType>(f);

      chamfer_distance_transform(f, nbh, out);

      return out;
    }


  }

}

#endif // ! MLN_TRANSFORM_CHAMFER_DISTANCE_TRANSFORM_HPP
