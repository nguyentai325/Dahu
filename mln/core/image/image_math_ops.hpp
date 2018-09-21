#ifndef MLN_CORE_IMAGE_IMAGE_MATH_OPS_HPP
# define MLN_CORE_IMAGE_IMAGE_MATH_OPS_HPP

# include <mln/core/image/image_expr.hpp>
# include <mln/core/math_ops.hpp>

namespace mln
{
  /**********************************/
  /***  Point-wise function       ***/
  /**********************************/

  MLN_GENERATE_CONST_UNARY_EXPR(abs, functional::abs_t);
  MLN_GENERATE_CONST_UNARY_EXPR(sqr, functional::sqr_t);
  MLN_GENERATE_CONST_UNARY_EXPR(pow, functional::pow_t);
  MLN_GENERATE_CONST_UNARY_EXPR(sqrt, functional::sqrt_t);
  MLN_GENERATE_CONST_UNARY_EXPR(cbrt, functional::cbrt_t);

  MLN_GENERATE_CONST_UNARY_EXPR(sum, functional::sum_t);
  MLN_GENERATE_CONST_UNARY_EXPR(prod, functional::prod_t);
  MLN_GENERATE_CONST_UNARY_EXPR(l0norm, functional::l0norm_t);
  MLN_GENERATE_CONST_UNARY_EXPR(l1norm, functional::l1norm_t);
  MLN_GENERATE_CONST_UNARY_EXPR(l2norm, functional::l2norm_t);
  MLN_GENERATE_CONST_UNARY_EXPR(linfnorm, functional::linfnorm_t);
}

#endif // ! MLN_CORE_IMAGE_IMAGE_MATH_OPS_HPP
