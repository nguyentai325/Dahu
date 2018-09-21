#ifndef MLN_MORPHO_TOS_IMPL_CTOS_SERIAL_HPP
# define MLN_MORPHO_TOS_IMPL_CTOS_SERIAL_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>
# include <mln/morpho/tos/pset.hpp>
# include <mln/morpho/tos/pset_priority.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      template <typename Compare>
      struct equiv;
    }


    namespace impl
    {

      namespace serial
      {


        template <typename I,
                  typename Neighborhood,
                  typename Compare,
                  typename Equiv,
                  bool use_priority = false>
        morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
        cToS(const Image<I>& ima,
             const Neighborhood& nbh,
             mln_point(I) pmin,
             const Compare& cmp,
             const Equiv& eq);


        /************************************/
        /**     Implementation             **/
        /************************************/


        namespace internal
        {
          template <typename I>
          inline
          typename I::size_type
          zfind_root(I& zpar, typename I::size_type x)
          {
            if (zpar[x] != x)
              zpar[x] = zfind_root(zpar, zpar[x]);
            return zpar[x];
          }

        }


        template <typename I,
                  typename Neighborhood,
                  typename Compare,
                  typename Equiv,
                  bool use_priority>
        morpho::component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
        cToS(const Image<I>& ima_,
             const Neighborhood& nbh,
             mln_point(I) pmin,
             const Compare& cmp,
             const Equiv& eq)
        {
          mln_entering("mln::morpho::impl::serial::ctos");
          using namespace mln::morpho::tos;

          typedef mln_value(I) V;
          typedef irange<V> R;
          typedef typename I::size_type size_type;

          static constexpr size_type UNPROCESSED = value_traits<size_type>::max();
          static constexpr size_type PROCESSED = 0;

          const I& ima = exact(ima_);


          // f: image of interval in Khalimsky space
          // K: image of value in Khalimsky that tells at which a level a point is inserted
          mln_ch_value(I, R) f = tos::internal::immerse(ima, cmp);
          mln_concrete(I) K(f, mln::init());


          component_tree<size_type, mln_ch_value(I, unsigned)> tree;
          auto& nodes = tree._get_data()->m_nodes;
          auto& S = tree._get_data()->m_S;
          auto& pmap = tree._get_data()->m_pmap;

          S.reserve(f.domain().size());
          pmap = imchvalue<unsigned>(f).init(UNPROCESSED);
          extension::fill(pmap, PROCESSED);

          typedef typename std::conditional<not use_priority,
                                            pset<mln_concrete(I), Compare>,
                                            pset_priority<mln_concrete(I), Compare> >::type pset_t;

        auto dindexes = wrt_delta_index(f, nbh.dpoints);

        // 1. Propagation
        {
          mln_entering("mln::morpho::impl::serial::ctos - propagation");

          pset_t W(K, cmp);
          size_type p = f.index_of_point(pmin * 2);
          W.insert(p);
          pmap[p] = PROCESSED;
          K[p] = f[p].lower;

          while (!W.empty())
            {
              p = W.has_next(p) ? W.pop_next(p) : W.pop_previous(p);
              V curlevel = K[p];
              S.push_back(p);

              mln_foreach (int k, dindexes)
                {
                  size_type q = p + k;
                  if (pmap[q] == UNPROCESSED)
                    {
                      if (cmp(f[q].upper, curlevel))
                        K[q] = f[q].upper;
                      else if (cmp(curlevel, f[q].lower))
                        K[q] = f[q].lower;
                      else
                        K[q] = curlevel;

                      pmap[q] = PROCESSED;
                      W.insert(q);
                    }
                }
            }
        }
        mln_exiting();
        // }

        // Warning parent & pmap refer to the same image
        // to avoid space waste.
        unsigned spos = S.size()-1;
        unsigned nnode = 0;
        mln_ch_value(I, unsigned)& parent = pmap;

        // 2nd step: union-find
        {
          mln_entering("mln::morpho::impl::serial::ctos - unionfind");
          mln_ch_value(I, unsigned) zpar;
          resize(zpar, K).init(UNPROCESSED);
          extension::fill(zpar, UNPROCESSED);

          auto is_face_2 = [](const point2d& p) { return p[0] % 2 == 0 and p[1] % 2 == 0; };

          for (int i = S.size()-1; i >= 0; --i)
            {
              size_type p = S[i];
              parent[p] = p;
              zpar[p] = p;
              nnode++;

              size_type rp = p;
              bool face2 = is_face_2(K.point_at_index(p));

              mln_foreach (int k, dindexes)
                {
                  size_type q = p + k;
                  if (zpar[q] != UNPROCESSED)
                    {
                      size_type r = internal::zfind_root(zpar, q);
                      if (r != rp) { // MERGE r and p
                        bool are_eq = eq(K[p], K[r]);
                        if (are_eq and !face2 and is_face_2(K.point_at_index(r)))
                          {
                            parent[rp] = r;
                            zpar[rp] = r;
                            face2 = true;
                            S[spos--] = rp;
                            rp = r;
                          }
                        else
                          {
                            parent[r] = rp;
                            zpar[r] = rp;
                            S[spos--] = r;
                          }
                        nnode -= are_eq;
                      }
                    }
                }
            }
          S[0] = internal::zfind_root(zpar, S[0]);
          mln_exiting();
        }


        // 4. Fully linked tree.
        //{
        mln_entering("mln::morpho::impl::serial::ctos - postprocessing");

        nodes.resize(nnode + 1);

        // Sentinel and Root
        nodes[0] = {0, 1, 0, 0, (unsigned)S.size()};
        nodes[1] = {0, 0, 0, 0, 0};
        pmap[S[0]] = 1;

        unsigned npos = 2;
        // 3rd step: canonicalization
        for (unsigned i = 1; i < S.size(); ++i)
          {
            size_type p = S[i];
            size_type q = parent[p];

            if (not eq(K[p], K[q])) // create a new node
              {
                pmap[p] = npos;
                unsigned par = pmap[q];
                unsigned next = nodes[par].m_next;
                nodes[npos].m_parent = par;

                // Insert node between par and next(par)
                nodes[npos].m_next = next;
                nodes[npos].m_prev = par;
                nodes[next].m_prev = npos;
                nodes[par].m_next = npos;

                // Update next sibling
                nodes[npos].m_next_sibling = next;

                nodes[npos].m_point_index = i;
                ++npos;
              }
            else
              {
                pmap[p] = pmap[q];
              }
          }
        nodes.resize(npos);
        mln_exiting();
        // }

        mln_postcondition(S.size() == K.domain().size());
        // All done !
        mln_exiting();
        return tree.get_subtree(1);
        }

      }

    }

  }

}


#endif // !MLN_MORPHO_TOS_IMPL_CTOS_SERIAL_HPP
