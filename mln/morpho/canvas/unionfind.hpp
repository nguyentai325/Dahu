#ifndef MLN_MORPHO_CANVAS_UNIONFIND_HPP
# define MLN_MORPHO_CANVAS_UNIONFIND_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/core/extension/extension.hpp>
# include <mln/core/algorithm/sort_sites.hpp>
# include <mln/core/dontcare.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace morpho
  {

    namespace canvas
    {

      /// \brief Define a default union-find visitor.
      struct default_unionfind_visitor
      {
        ///  \brief Called on every point \p during the make-set step.
        void on_make_set(dontcare_t p) { (void) p;}

        /// \brief Called during the merge-step.
        ///
        /// Called when a component rooted in \p p merges with a
        /// component rooted in \p q, \p q becomes the new root.
        void on_union(dontcare_t p, dontcare_t q) { (void) p; (void) q; };

        /// \brief Called at the end of the algorithm for each point
        ///
        /// \p p is any point in the domain, \p q is the root of the
        /// component \p p belongs to.
        void on_finish(dontcare_t p, dontcare_t q) { (void) p; (void)q;};
      };

      template <class I, class N, class StopCriterion,
                class Compare = std::less<mln_value(I)>,
                class uf_visitor = default_unionfind_visitor>
      mln_ch_value(I, mln_point(I))
      unionfind(const Image<I>& input, const Neighborhood<N>& nbh,
                StopCriterion stop,
                Compare cmp = Compare(),
                uf_visitor viz = uf_visitor());

      /*********************/
      /** Implementation  **/
      /*********************/

      namespace internal
      {
        enum e_unionfind_status {
          FAIL = 0,
          PASS = 1,
          NONE = (unsigned char)-1
        };

        template <class I>
        mln_point(I)
        zfindroot(I& par, const mln_point(I)& p)
        {
          mln_point(I) q = par(p);
          if (q != p)
            return (par(p) = zfindroot(par, q));
          else
            return p;
        }


        template <class I, class N, class uf_visitor, class Compare, class StopCriterion, class J>
        void
        unionfind_impl(const I& input, const N& nbh, StopCriterion term, Compare cmp, uf_visitor viz,
                         mln_ch_value(I, mln_point(I))& par, J&& status)
        {
          auto S = sort_sites(input, cmp);

          mln_point(I) p;
          mln_iter(nit, nbh(p));
          mln_foreach (p, S)
            {
              // make-set
              par(p) = p;
              viz.on_make_set(p);
              status(p) = term(p);

              mln_foreach(mln_point(I) n, nit)
                if (status.at(n) != NONE)
                  {
                    mln_point(I) r = zfindroot(par, n);
                    // merge step
                    if (p != r)
                      {
                        if (status(r) == PASS) // r pass already, no need to merge.
                          {
                            status(p) = PASS;
                          }
                        else
                          {
                            viz.on_union(r, p);
                            par(r) = p;
                            status(p) |= status(r);
                          }
                      }
                  }

              status(p) = status(p) || term(p);
            }
          // Set the root status to PASS
          status(S.back()) = PASS;

          // Reverse, canonize and call on_finish
          mln_reverse_foreach(p, S)
            {
              mln_point(I) q = par(p);
              par(p) = par(q);
              assert(status(par(p)));

              viz.on_finish(p, par(p));
            }

        }


        //
        // The facade chooses the best way to proceed the union-find
        // 1. If the input image is able to hold the neighborhood
        //    nothing special occurs.
        // 2. Otherwise, a fake extension is added to the status image
        template <class I, class N, class uf_visitor, class Compare, class StopCriterion>
        mln_ch_value(I, mln_point(I))
          unionfind_facade(const I& input, const N& nbh, StopCriterion term, Compare cmp, uf_visitor viz)
        {
          mln_ch_value(I, mln_point(I)) par = imchvalue<mln_point(I)>(input);
          mln_ch_value(I, unsigned char) status = imchvalue<unsigned char>(input).init(NONE);

          if (not extension::need_adjust(status, nbh))
            unionfind_impl(input, nbh, term, cmp, viz, par, status);
          else
            {
              mln::trace::warn("Slow version because input image extension is not wide enough.");
              unionfind_impl(input, nbh, term, cmp, viz, par,
                             extension::add_value_extension(status, NONE));
            }

          return par;
        }
      }

      template <class I, class N, class StopCriterion, class Compare, class uf_visitor>
      mln_ch_value(I, mln_point(I))
      unionfind(const Image<I>& input,
                const Neighborhood<N>& nbh,
                StopCriterion stop,
                Compare cmp,
                uf_visitor viz)
      {
        return internal::unionfind_facade(exact(input), exact(nbh), stop, cmp, viz);
      }


      /*********************************************/
      /***       Visitors                       ****/
      /*********************************************/

    }

  }

}

#endif // ! MLN_MORPHO_CANVAS_UNIONFIND_HPP
