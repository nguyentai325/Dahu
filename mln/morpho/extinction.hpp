#ifndef MLN_MORPHO_EXTINCTION_HPP
# define MLN_MORPHO_EXTINCTION_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/morpho/canvas/unionfind.hpp>
# include <mln/core/always.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace morpho
  {


    /// \brief Compute the extinction value of each minima
    ///
    /// If a node is not a local minima, its extinction value is null.
    /// Otherwise, it is the value of the attribute for which the minima
    /// is extincted. If no attribute is given, this is the dynamic.
    template <class I, class N, class Compare = std::less<mln_value(I)> >
    mln_concrete(I)
    extinction(const Image<I>& ima,
               const Neighborhood<N>& nbh,
               const Compare& cmp = Compare());

    /***********************/
    /*** Implementation  ***/
    /***********************/

    namespace internal
    {

      template <class I, class Compare>
      struct extinction_dynamic_ufind_visitor
      {

        void on_make_set(const mln_point(I)& p)
        {
          amin(p) = p;
        }

        void on_union(const mln_point(I)& p,
                      const mln_point(I)& q)
        {
          mln_point(I) mp = amin(p);
          mln_point(I) mq = amin(q);

          if (cmp(f(mp), f(mq))) {      // mp is the minimum
            amin(q) = mp;
            extinction(mq) = abs(f(q) - f(mq));
          } else {                      // mq is the minimum
            amin(q) = mq;
            extinction(mp) = abs(f(q) - f(mp));
          }
        }

        void on_finish(dontcare_t, dontcare_t)
        {
        }

        const I&                        f;
        mln_ch_value(I, mln_point(I))&  amin;
        mln_concrete(I)&                extinction;
        Compare                         cmp;
      };

    }

    template <class I, class N, class Compare>
    mln_concrete(I)
    extinction(const Image<I>& ima_,
               const Neighborhood<N>& nbh_,
               const Compare& cmp)
    {
      mln_entering("mln::morpho::extinction");

      const I& ima = exact(ima_);
      const N& nbh = exact(nbh_);

      mln_concrete(I) extinction = imconcretize(ima);
      mln_ch_value(I, mln_point(I)) amin = imchvalue<mln_point(I)>(ima);


      internal::extinction_dynamic_ufind_visitor<I, Compare>
        viz = {ima, amin, extinction, cmp};

      auto par = canvas::unionfind(ima, nbh, no, cmp, viz);
      // Set root extinction value
      {
        mln_viter(v, par);
        v.init();
        mln_point(I) root = *v;
        extinction(amin(root)) = abs(ima(root) - ima(amin(root)));
      }

      // Propagate extinction value to minima flatzones ?
      // Is that necessary ?

      mln_exiting();
      return extinction;
    }


  }

}

#endif // ! MLN_MORPHO_EXTINCTION_HPP
