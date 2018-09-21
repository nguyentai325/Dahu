#ifndef MLN_CORE_VEC_VEC_MATH_OPS_HPP
# define MLN_CORE_VEC_VEC_MATH_OPS_HPP

# include <mln/core/vec.hpp>
# include <mln/core/math_ops.hpp>


/**
* \file Define vectorial overloads for fundamental mathematical and
* statistical operators
*
*/

# define MLN_VEC_PROMOTE_FUN(T, dim, tag, fun)				\
  internal::vec_base< decltype(fun(std::declval<T>())), dim, tag>

# define MLN_PROMOTE(T)					\
  decltype(std::declval<T>() + std::declval<T>())

# define MLN_PROMOTE_FUN(T, fun)		\
  decltype(fun(std::declval<T>()))

# define MLN_PROMOTE_FUNCOMP(T, f, g)		\
  decltype(f(g(std::declval<T>())))

namespace mln
{

  /* Element wise operators */
  template <class U, class V, unsigned dim, class tag>
  internal::vec_base< decltype( pow(std::declval<U>(), std::declval<V>()) ),
                      dim, tag>
  pow(const internal::vec_base<U, dim, tag>& x, V exp);

  template <class T, unsigned dim, class tag>
  MLN_VEC_PROMOTE_FUN(T, dim, tag, sqr)
    sqr(const internal::vec_base<T, dim, tag>& x);

  template <class T, unsigned dim, class tag>
  MLN_VEC_PROMOTE_FUN(T, dim, tag, abs)
    abs(const internal::vec_base<T, dim, tag>& x);

  template <class T, unsigned dim, class tag>
  MLN_VEC_PROMOTE_FUN(T, dim, tag, sqrt)
    sqrt(const internal::vec_base<T, dim, tag>& x);

  template <class T, unsigned dim, class tag>
  MLN_VEC_PROMOTE_FUN(T, dim, tag, cqrt)
    cbrt(const internal::vec_base<T, dim, tag>& x);


  /* Reduction operators */

  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE(T)
  sum(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE(T)
  prod(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  T
  maximum(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  T
  minimum(const internal::vec_base<T, dim, tag>& x);


  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUN(T, abs)
    l0norm(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUN(T, abs)
    linfnorm(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUN(T, abs)
    l1norm(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUNCOMP(T, sqrt, sqr)
    l2norm(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUN(T, sqr)
    l2norm_sqr(const internal::vec_base<T, dim, tag>& x);

  template <unsigned p, typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUNCOMP(T, pow, abs)
    lpnorm(const internal::vec_base<T, dim, tag>& x);

  template <typename T, unsigned dim, class tag>
  inline
  auto
  l2dist(const internal::vec_base<T, dim, tag>& x,
         const internal::vec_base<T, dim, tag>& y)
    -> decltype(l2norm(x-y));


  /*****************************/
  /**   Implementation       ***/
  /*****************************/

# define MLN_GEN_CODE(FUN)			\
  template <class T, unsigned dim, class tag>	\
  inline					\
  MLN_VEC_PROMOTE_FUN(T, dim, tag, FUN)		\
    FUN(const internal::vec_base<T, dim, tag>& x)		\
  {						\
    MLN_VEC_PROMOTE_FUN(T, dim, tag, FUN) res;	\
    for (unsigned i = 0; i < dim; ++i)		\
      res[i] = FUN(x[i]);			\
    return res;					\
  }

#define MLN_GEN_BINARY_CODE(FUN, FUNBASE, EXPR)                  \
  template <class T, unsigned dim, class tag>                    \
  inline                                                         \
  auto                                                           \
  FUN(const internal::vec_base<T, dim, tag>& x,                  \
      const internal::vec_base<T, dim, tag>& y)                  \
    -> decltype(FUNBASE(EXPR))                                   \
  {                                                              \
    return FUNBASE(EXPR);                                        \
  }



  MLN_GEN_CODE(sqr);
  MLN_GEN_CODE(sqrt);
  MLN_GEN_CODE(abs);
  MLN_GEN_CODE(cbrt);

  MLN_GEN_BINARY_CODE(l2dist, l2norm, x-y);


  template <class U, class V, unsigned dim, class tag>
  internal::vec_base< decltype( pow(std::declval<U>(), std::declval<V>()) ),
                      dim, tag>
  pow(const internal::vec_base<U, dim, tag>& x, V exp)
  {
    typedef decltype(pow(std::declval<U>(), std::declval<V>())) R;
    internal::vec_base<R, dim, tag> res;
    for (unsigned i = 0; i < dim; ++i)
      res[i] = pow(x[i], exp);
    return res;
  }


  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE(T)
  sum(const internal::vec_base<T, dim, tag>& x)
  {
    MLN_PROMOTE(T) res = x[0];
    for (unsigned i = 1; i < dim; ++i)
      res += x[i];
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE(T)
  prod(const internal::vec_base<T, dim, tag>& x)
  {
    MLN_PROMOTE(T) res = x[0];
    for (unsigned i = 1; i < dim; ++i)
      res *= x[i];
    return res;
  }


  template <typename T, unsigned dim, class tag>
  inline
  T
  maximum(const internal::vec_base<T, dim, tag>& x)
  {
    T res = x[0];
    for (unsigned i = 1; i < dim; ++i)
      if (res < x[i])
	res = x[i];
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  T
  minimum(const internal::vec_base<T, dim, tag>& x)
  {
    T res = x[0];
    for (unsigned i = 1; i < dim; ++i)
      if (x[i] < res)
	res = x[i];
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUN(T, abs)
    l0norm(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUN(T, abs) U;
    U res = abs(x[0]);
    for (unsigned i = 1; i < dim; ++i)
      {
	U v = abs(x[i]);
	if (v < res)
	  res = v;
      }
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUN(T, abs)
    linfnorm(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUN(T, abs) U;
    U res = abs(x[0]);
    for (unsigned i = 1; i < dim; ++i)
      {
	U v = abs(x[i]);
	if (res < v)
	  res = v;
      }
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUN(T, abs)
    l1norm(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUN(T, abs) U;
    U res = abs(x[0]);
    for (unsigned i = 1; i < dim; ++i)
      res += abs(x[i]);
    return res;
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUNCOMP(T, sqrt, sqr)
    l2norm(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUNCOMP(T, sqrt, sqr) U;
    U res = sqr(x[0]);
    for (unsigned i = 1; i < dim; ++i)
      res += sqr(x[i]);
    return sqrt(res);
  }

  template <typename T, unsigned dim, class tag>
  inline
  MLN_PROMOTE_FUN(T, sqr)
    l2norm_sqr(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUNCOMP(T, sqrt, sqr) U;
    U res = sqr(x[0]);
    for (unsigned i = 1; i < dim; ++i)
      res += sqr(x[i]);
    return res;
  }


  template <unsigned p, typename T, unsigned dim, class tag>
  MLN_PROMOTE_FUNCOMP(T, pow, abs)
    lpnorm(const internal::vec_base<T, dim, tag>& x)
  {
    typedef MLN_PROMOTE_FUNCOMP(T, pow, abs) U;
    U res = pow(abs(x[0]), p);
    for (unsigned i = 1; i < dim; ++i)
      res += pow(abs(x[i]), p);
    return pow(res, 1/p);
  }

}

# undef MLN_GEN_CODE
# undef MLN_GEN_BINARY_CODE
# undef MLN_VEC_PROMOTE_FUN
# undef MLN_VEC_PROMOTE
# undef MLN_PROMOTE_FUNCOMP

#endif // ! MLN_CORE_VEC_VEC_MATH_OPS_HPP
