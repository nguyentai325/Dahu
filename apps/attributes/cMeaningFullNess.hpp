#ifndef APPS_ATTRIBUTE_CMEANINGFULLNESS_HPP
# define APPS_ATTRIBUTE_CMEANINGFULLNESS_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/neighborhood/dyn_neighborhood.hpp>
# include <mln/core/trace.hpp>
# include <mln/accu/accumulator.hpp>
# include <mln/accu/accumulators/variance.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>
# include <mln/morpho/component_tree/compute_depth.hpp>
# include <apps/tos/topology.hpp>
# include <apps/tos/croutines.hpp>
# include <apps/g2/accu/lca.hpp>
# include "curvature.hpp"

namespace mln
{

  /// \brief compute the external regional energy
  ///
  template <class P, class AMap, class I>
  property_map<morpho::component_tree<P, AMap>, float>
  compute_meaningfullness_external_energy(const morpho::component_tree<P, AMap>& tree,
                                          const Image<I>& vmap,
                                          int eps);

  /// \brief generic function to compute a contextual energy
  /// This function supposes a tree computed with 0 and 1-faces. The accumulation
  /// is performed on 2F only.
  /// Considering a ball Β of radius ε, the accumulator is computed for every shape Γ on:
  ///  * Ext(Γ) = δ(Γ) \ Γ
  ///  * Int(Γ) = Γ \ ε(Γ)
  ///  * Reg(Γ) = Ext(Γ) ⋃ Int(Γ)
  ///
  /// \param tree The component tree.
  /// \param valuemap Image with values to accumulate (only 2F are used).
  /// \param accu The accumulator (-like) to compute.
  /// \param eps The radius of the ball used as the SE for shape dilation and erosion.
  ///
  /// \return An attribute map for each node with the 3 values of the
  ///         accumulation on thes internal/external/whole regions.
  template <class P, class I, class AccuLike>
  property_map<morpho::component_tree<P, image2d<P> >,
               vec<typename accu::result_of<AccuLike, mln_value(I)>::type, 3> >
  compute_regional_energy(const morpho::component_tree<P, image2d<P> >& tree,
                          const Image<I>& valuemap,
                          const AccumulatorLike<AccuLike>& accu,
                          int eps);


  /*********************/
  /** Implementation ***/
  /*********************/

  template <class P, class I, class AccuLike>
  property_map<morpho::component_tree<P, image2d<P> >,
               vec<typename accu::result_of<AccuLike, mln_value(I)>::type, 3> >
  compute_regional_energy(const morpho::component_tree<P, image2d<P> >& tree,
                          const Image<I>& valuemap,
                          const AccumulatorLike<AccuLike>& accu_,
                          int eps)
  {
    mln_entering("mln::compute_regional_energy");

    typedef morpho::component_tree<P, image2d<P> > tree_t;
    typedef typename tree_t::node_type node_type;
    typedef typename tree_t::vertex_id_t vertex_id_t;
    typedef typename accu::result_of<AccuLike, mln_value(I)>::type R;

    const I& vmap = exact(valuemap);
    auto acc = accu::make_accumulator(exact(accu_), mln_value(I) ());

    typedef decltype(acc) Accu;

    // Check that domain matches.
    mln_precondition(tree._get_data()->m_pmap.domain() == vmap.domain());


    acc.init();
    property_map<tree_t, Accu> external(tree, acc);
    property_map<tree_t, Accu> internal(tree, acc);

    // Create the se (ball only on 2F)
    std::vector<point2d>      dpoints;
    for (int i = -eps; i < eps+1; ++i)
      for (int j = -eps; j < eps+1; ++j)
        if ( (i != 0 or j != 0) and (i*i + j*j) <= (eps*eps) )
          dpoints.emplace_back(2*i,2*j);

    dyn_neighborhood<std::vector<point2d>,
                     dynamic_neighborhood_tag>  ball(dpoints);


    // Compute the depth image.
    property_map<tree_t, unsigned> depth = morpho::compute_depth(tree);

    // Compute internal region attribute
    {
      mln_pixter(px, vmap);
      mln_iter(nx, ball(px));
      accu::least_common_ancestor<P, image2d<P>>  lca(tree, depth);

      mln_forall(px)
      {
        point2d p = px->point();
        if (K1::is_face_2(p))
          {
            lca.init();
            lca.take(px->index());

            mln_forall(nx)
              if (vmap.domain().has(nx->point()))
                lca.take(nx->index());
              else
                break; // border pixel LCA is above the root

            node_type x = tree.get_node_at(px->index());
            node_type y = nx.finished() ? lca.to_result() : tree.nend();
            for (;x != y; x = x.parent())
              internal[x].take(px->val());
          }
      }
    }

    // Compute external region attribute
    {
      property_map<tree_t, bool> dejavu(tree, false);
      std::vector< std::pair<vertex_id_t, vertex_id_t> > branches;

      mln_pixter(px, vmap);
      mln_iter(nx, ball(px));
      mln_forall(px)
      {
        point2d p = px->point();
        if (K1::is_face_2(p))
          {
            node_type k = tree.get_node_at(px->index());
            branches.clear();

            mln_forall(nx)
              if (vmap.domain().has(nx->point()))
                {
                  node_type x = tree.get_node_at(nx->index());
                  node_type y = lca(tree, depth, k, x);
                  branches.emplace_back(x.id(), y.id());
                  for (node_type n = x; n != y; n = n.parent())
                    if (not dejavu[n]) {
                      external[n].take(px->val());
                      dejavu[n] = true;
                    }
                }

            // undo dejavu
            {
              vertex_id_t x,y;
              mln_foreach(std::tie(x,y), branches)
                for (node_type n = tree.get_node(x); n.id() != y; n = n.parent())
                  dejavu[n] = false;
            }
          }
      }
    }

    // Return the results.
    property_map<tree_t, vec<R, 3> > res(tree);
    {
      mln_foreach(auto x, tree.nodes())
        {
          res[x][0] = internal[x].to_result();
          res[x][1] = external[x].to_result();
          internal[x].take(external[x]);
          res[x][2] = internal[x].to_result();
        }
    }

    mln_exiting();
    return res;
  }


  template <class P, class AMap, class I>
  property_map<morpho::component_tree<P, AMap>, float>
  compute_meaningfullness_external_energy(const morpho::component_tree<P, AMap>& tree,
                                          const Image<I>& vmap,
                                          int eps)
  {
    auto attrmap = compute_regional_energy(tree, vmap,
                                           accu::features::variance<double> (),
                                           eps);

    property_map<morpho::component_tree<P, AMap>, float> energy(tree);
    mln_foreach(auto node, tree.nodes())
      energy[node] = (attrmap[node][0] + attrmap[node][1]) / attrmap[node][2];

    return energy;
  }

}

#endif // ! APPS_ATTRIBUTE_CMEANINGFULLNESS_HPP
