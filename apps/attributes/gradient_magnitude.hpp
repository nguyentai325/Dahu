#ifndef APPS_ATTRIBUTE_GRADIENT_MAGNITUDE_HPP
# define APPS_ATTRIBUTE_GRADIENT_MAGNITUDE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/math_ops.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <apps/tos/croutines.hpp>

namespace mln
{

  template <class P, class AMap, class V>
  property_map<morpho::component_tree<P, AMap>, float>
  compute_gradient_magnitude(const morpho::component_tree<P, AMap>& tree,
                             const image2d<V>& f)
  {
    mln_entering("mln::compute_gradient_magnitude");

    image2d<float> grad;
    resize(grad, f).init(1.0);

    box2d dom = f.domain();
    sbox2d d2f = {dom.pmin - point2d{2,2}, dom.pmax, point2d{2,2}};

    float lmax = l1norm(value_traits<V>::sup() - value_traits<V>::inf());

    mln_foreach(point2d p, d2f)
      {
        grad.at(p + point2d{0,1}) = 1 - l1norm(f.at(p) - f.at(p + point2d{0,2})) / lmax;
        grad.at(p + point2d{1,0}) = 1 - l1norm(f.at(p) - f.at(p + point2d{2,0})) / lmax;
        grad.at(p + point2d{1,1}) = 1 - l1norm(f.at(p) - f.at(p + point2d{2,2})) / lmax;
      }

    auto res = compute_attribute_on_contour(tree, grad, accu::features::mean<> ());
    res[tree.get_root()] = 1.0;
    mln_exiting();
    return res;
  }

}


#endif // ! APPS_ATTRIBUTE_GRADIENT_MAGNITUDE_HPP
