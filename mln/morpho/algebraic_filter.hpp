#ifndef MLN_MORPHO_ALGEBRAIC_FILTER_HPP
# define MLN_MORPHO_ALGEBRAIC_FILTER_HPP

# include <mln/accu/accumulators/count.hpp>
# include <mln/morpho/canvas/unionfind.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace morpho
  {

    /// \brief Compute the area algebraic closing of an image.
    /// \param ima The input image
    /// \param nbh The neighborhood
    /// \param area The grain size
    /// \param cmp A strict total ordering on values
    template <class I, class N, class Compare = std::less<mln_value(I)> >
    mln_concrete(I)
    area_closing(const Image<I>& ima,
                 const Neighborhood<N>& nbh,
                 unsigned area,
                 Compare cmp = Compare());

    /// \brief Compute the area algebraic opening of an image.
    /// \param ima The input image
    /// \param nbh The neighborhood
    /// \param area The grain size
    template <class I, class N>
    mln_concrete(I)
    area_opening(const Image<I>& ima,
                 const Neighborhood<N>& nbh,
                 unsigned area);



    /******************************/
    /*** Implementation          **/
    /******************************/

    namespace internal
    {
      template <class I, class Accu>
      struct attribute_and_reconstruct_ufind_visitor
      {
        void on_make_set(const mln_point(I)& p)
        {
          m_acc_img(p).take(p);
        }

        void on_union(const mln_point(I)& p, const mln_point(I)& q)
        {
          m_acc_img(q).take(m_acc_img(p));
        }

        void on_finish(const mln_point(I)& p, const mln_point(I)& q)
        {
          m_rec(p) = m_rec(q);
        }

        mln_ch_value(I, Accu)&    m_acc_img;
        mln_concrete(I)&          m_rec;
      };

    }

    template <class I, class N, class Compare>
    mln_concrete(I)
    area_closing(const Image<I>& ima_,
                 const Neighborhood<N>& nbh,
                 unsigned area,
                 Compare cmp)
    {
      mln_entering("mln::morpho::area_closing");

      const I& ima = exact(ima_);

      typedef mln::accu::accumulators::count<typename I::size_type> ACCU;

      mln_ch_value(I, ACCU) accu_img = imchvalue<ACCU>(ima);
      mln_concrete(I) out = clone(ima);

      internal::attribute_and_reconstruct_ufind_visitor<I, ACCU> viz {accu_img, out};
      mln::morpho::canvas::unionfind(ima, nbh,
                                     [&accu_img, area] (const mln_point(I)& p) {
                                       return accu_img(p).to_result() >= area; },
                                     cmp, viz);
      mln_exiting();
      return out;
    }

    template <class I, class N>
    mln_concrete(I)
    area_opening(const Image<I>& ima_,
                 const Neighborhood<N>& nbh,
                 unsigned area)
    {
      mln_entering("mln::morpho::area_opening");
      mln_concrete(I) out = area_closing(exact(ima_), exact(nbh), area,
                                         std::greater<mln_value(I)> ());
      mln_exiting();
      return out;
    }

  }

}

#endif // ! MLN_MORPHO_ALGEBRAIC_FILTER_HPP
