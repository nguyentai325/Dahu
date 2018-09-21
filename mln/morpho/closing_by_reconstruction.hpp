#ifndef MLN_MORPHO_CLOSING_BY_RECONSTRUCTION_HPP
# define MLN_MORPHO_CLOSING_BY_RECONSTRUCTION_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/morpho/canvas/unionfind.hpp>
# include <mln/accu/accumulators/min.hpp>

namespace mln
{

  namespace morpho
  {

    template <class I, class J, class N, class Compare = std::less<mln_value(I)> >
    mln_concrete(I)
    closing_by_reconstruction(const Image<I>& f,
                              const Image<J>& markers,
                              const Neighborhood<N>& nbh,
                              Compare cmp = Compare());


    /*************************/
    /** Implementation     ***/
    /*************************/

    namespace internal
    {

      template <class I, class J, class Compare>
      struct closing_by_rec_ufind_visitor
      {
        typedef mln::accu::accumulators::min<mln_value(J), Compare> accu_t;

        closing_by_rec_ufind_visitor(const J& markers, mln_ch_value(I, accu_t)& accus, mln_concrete(I)& output)
          : m_markers(markers),
            m_accu (accus),
            m_out(output)
        {
        }


        void on_make_set(const mln_point(I)& p)
        {
          m_accu(p).take(m_markers(p));
        }

        void on_union(const mln_point(I)& p, const mln_point(I)& q)
        {
          m_accu(q).take(m_accu(p));
        }

        void on_finish(const mln_point(I)& p, mln_point(I)& q)
        {
          m_out(p) = m_out(q);
        }

        const J&                 m_markers;
        mln_ch_value(J, accu_t)& m_accu;
        mln_concrete(I)&         m_out;
      };

    }

    template <class I, class J, class N, class Compare>
    mln_concrete(I)
    closing_by_reconstruction(const Image<I>& f_,
                              const Image<J>& markers_,
                              const Neighborhood<N>& nbh,
                              Compare cmp)
    {
      static_assert( std::is_convertible<mln_value(J), mln_value(I)>::value,
                     "Marker image value type must be convertible to f's value type.");

      static_assert( std::is_same<mln_point(J), mln_point(I)>::value,
                     "Images f and marker must have the same point type." );

      mln_entering("mln::morpho::closing_by_reconstruction");
      typedef mln::accu::accumulators::min<mln_value(J), Compare> accu_t;
      const I&                  f = exact(f_);
      const J&                  markers = exact(markers_);

      // assert that f <= markers
      mln_precondition(all(imtransform( imzip(f,markers), [cmp](std::tuple<mln_value(I), mln_value(J)> v) {
              return not cmp(std::get<1>(v), std::get<0>(v));
            })));

      mln_concrete(I)           out = clone(f);
      mln_ch_value(I, accu_t)   accus;

      resize(accus, f).init( accu_t(cmp) );
      internal::closing_by_rec_ufind_visitor<I, J, Compare> viz(markers, accus, out);

      auto criterion = [&accus, &f, cmp] (const mln_point(I)& p) {
        return not cmp(f(p), accus(p).to_result());
      };

      morpho::canvas::unionfind(f, nbh, criterion, cmp, viz);
      mln_exiting();
      return out;
    }

  }

}

#endif // ! MLN_MORPHO_CLOSING_BY_RECONSTRUCTION_HPP
