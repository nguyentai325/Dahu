#ifndef APPS_ATTRIBUTE_CMSER_HPP
# define APPS_ATTRIBUTE_CMSER_HPP

# include <mln/core/trace.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>

enum eMSER_attribute {
  MSER_DIFF,   // Val > 0,      δΓ = A(Γ⁺) - A(Γ⁻)
  MSER_RATIO,  // 0 < Val < 1   δΓ = A(Γ⁻) / A(Γ⁺)
  MSER_NORM    // 0 < Val < 1   δΓ = (A(Γ⁺) - A(Γ⁻)) / A(Γ)
};

enum eMSER_accum_type {
  MSER_ABSOLUTE,  // Γ⁺ = ⋁ { Γ ⊂ Γ', |V(Γ') - V(Γ)| ≥ ε }
  MSER_SUM        // Γ⁺ = ⋁ { Γ ⊂ Γ', ∑ |V(par(u)) - V(u)| > ε, u ∈ [Γ → Γ'] }
};

namespace mln
{

  template <class P, class Amap,
            class ValuePropertyMap,
            class AreaPropertyMap,
            class Distance>
  property_map<morpho::component_tree<P, Amap>, float>
  compute_MSER(const morpho::component_tree<P, Amap>& tree,
               const ValuePropertyMap& vmap,
               const AreaPropertyMap& amap,
               float eps,
               eMSER_attribute amser = MSER_DIFF,
               eMSER_accum_type fsum = MSER_ABSOLUTE,
               Distance dist = Distance());


  template <class P, class Amap,
            class ValuePropertyMap,
            class AreaPropertyMap>
  property_map<morpho::component_tree<P, Amap>, float>
  compute_MSER(const morpho::component_tree<P, Amap>& tree,
               const ValuePropertyMap& vmap,
               const AreaPropertyMap& amap,
               float eps,
               eMSER_attribute amser = MSER_DIFF,
               eMSER_accum_type fsum = MSER_ABSOLUTE);


  /**********************/
  /**  Implementation  **/
  /**********************/

  template <class P, class Amap,
            class ValuePropertyMap,
            class AreaPropertyMap,
            class Distance>
  property_map<morpho::component_tree<P, Amap>, float>
  compute_MSER(const morpho::component_tree<P, Amap>& tree,
               const ValuePropertyMap& vmap,
               const AreaPropertyMap& amap,
               float eps,
               eMSER_attribute amser,
               eMSER_accum_type fsum,
               Distance dist)
  {
    mln_entering("mln::compute_MSER");

    typedef morpho::component_tree<P, Amap> tree_t;
    typedef typename tree_t::node_type node_t;

    property_map<tree_t, unsigned> aplus(tree);
    property_map<tree_t, unsigned> aminus(tree, 0);

    if (fsum == MSER_ABSOLUTE)
      {
        mln_reverse_foreach(auto n, tree.nodes())
          {
            auto x = n;
            while (x.get_parent_id() != tree.npos() and dist(vmap[n], vmap[x]) < eps)
              x = x.parent();

            aminus[x] = std::max<unsigned>(aminus[x], amap[n]);
            aplus[n] = amap[x];
          }
      }
    else
      {
        mln_reverse_foreach(auto n, tree.nodes())
          {
            auto x = n;
            float d = 0;
            while (x.id() != x.get_parent_id() and d < eps)
              {
                d += dist(vmap[x], vmap[x.parent()]);
                x = x.parent();
              }
            aminus[x] = std::max<unsigned>(aminus[x], amap[n]);
            aplus[n] = amap[x];
          }
      }
    aplus[tree.get_root_id()] = amap[tree.get_root()];

    std::function<float(const node_t&)> f;
    switch (amser) {
      case MSER_DIFF:
        f = [&](const node_t& x) -> float { return aplus[x] - aminus[x]; };
        break;
      case MSER_RATIO:
        f = [&](const node_t& x) -> float { return (float)amap[x] / aplus[x]; };
        break;
      case MSER_NORM:
        f = [&](const node_t& x) -> float { return (float)(aplus[x] - aminus[x]) / amap[x]; };
        break;
    };

    property_map<tree_t, float> out(tree);
    mln_foreach(auto x, tree.nodes())
      out[x] = f(x);

    mln_exiting();
    return out;
  }

  template <class P, class Amap,
            class ValuePropertyMap,
            class AreaPropertyMap>
  property_map<morpho::component_tree<P, Amap>, float>
  compute_MSER(const morpho::component_tree<P, Amap>& tree,
               const ValuePropertyMap& vmap,
               const AreaPropertyMap& amap,
               float eps,
               eMSER_attribute amser,
               eMSER_accum_type fsum)
  {
    typedef typename ValuePropertyMap::value_type V;
    auto dist = [](V x, V y) -> float { return l2norm(x - y); };

    return compute_MSER(tree, vmap, amap, eps, amser, fsum, dist);
  }

}

#endif // ! APPS_ATTRIBUTE_CMSER_HPP
