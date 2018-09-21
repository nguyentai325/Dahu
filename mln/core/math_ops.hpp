#ifndef MLN_CORE_MATH_OPS_HPP
# define MLN_CORE_MATH_OPS_HPP

/**
* \file
* \brief Define fundamental mathematical and statistical operators
*
*/

# include <utility>
# include <cmath>

namespace mln
{
  /* Element wise operators */
  int			sqr(int x);
  long			sqr(long n);
  long long		sqr(long long n);
  unsigned int		sqr(unsigned int x);
  unsigned long		sqr(unsigned long n);
  unsigned long long	sqr(unsigned long long n);
  float			sqr(float x);
  double		sqr(double x);
  long double		sqr(long double x);

  int			abs(int x);
  long			abs(long n);
  long long		abs(long long n);
  unsigned int		abs(unsigned int x);
  unsigned long		abs(unsigned long n);
  unsigned long long	abs(unsigned long long n);
  float			abs(float x);
  double		abs(double x);
  long double		abs(long double x);


  using std::sqrt;
  using std::cbrt;
  using std::pow;

  namespace functional
  {
    template <typename T = void> struct sqrt_t;
    template <typename T = void> struct cbrt_t;
    template <typename T = void> struct pow_t;
    template <typename T = void> struct sqr_t;
    template <typename T = void> struct abs_t;
  };


  /* Reduction operators */
  int			sum(int x);
  long			sum(long n);
  long long		sum(long long n);
  unsigned int		sum(unsigned int x);
  unsigned long		sum(unsigned long n);
  unsigned long long	sum(unsigned long long n);
  float			sum(float x);
  double		sum(double x);
  long double		sum(long double x);

  int			prod(int x);
  long			prod(long n);
  long long		prod(long long n);
  unsigned int		prod(unsigned int x);
  unsigned long		prod(unsigned long n);
  unsigned long long	prod(unsigned long long n);
  float			prod(float x);
  double		prod(double x);
  long double		prod(long double x);

  int			minimum(int x);
  long			minimum(long n);
  long long		minimum(long long n);
  unsigned int		minimum(unsigned int x);
  unsigned long		minimum(unsigned long n);
  unsigned long long	minimum(unsigned long long n);
  float			minimum(float x);
  double		minimum(double x);
  long double		minimum(long double x);

  int			maximum(int x);
  long			maximum(long n);
  long long		maximum(long long n);
  unsigned int		maximum(unsigned int x);
  unsigned long		maximum(unsigned long n);
  unsigned long long	maximum(unsigned long long n);
  float			maximum(float x);
  double		maximum(double x);
  long double		maximum(long double x);

  int			l0norm(int x);
  long			l0norm(long n);
  long long		l0norm(long long n);
  unsigned int		l0norm(unsigned int x);
  unsigned long		l0norm(unsigned long n);
  unsigned long long	l0norm(unsigned long long n);
  float			l0norm(float x);
  double		l0norm(double x);
  long double		l0norm(long double x);

  int			l1norm(int x);
  long			l1norm(long n);
  long long		l1norm(long long n);
  unsigned int		l1norm(unsigned int x);
  unsigned long		l1norm(unsigned long n);
  unsigned long long	l1norm(unsigned long long n);
  float			l1norm(float x);
  double		l1norm(double x);
  long double		l1norm(long double x);

  int			l2norm(int x);
  long			l2norm(long n);
  long long		l2norm(long long n);
  unsigned int		l2norm(unsigned int x);
  unsigned long		l2norm(unsigned long n);
  unsigned long long	l2norm(unsigned long long n);
  float			l2norm(float x);
  double		l2norm(double x);
  long double		l2norm(long double x);

  int			l2dist(int x, int y);
  long			l2dist(long x, long y);
  long long		l2dist(long long x, long long y);
  unsigned int		l2dist(unsigned int x, unsigned int y);
  unsigned long		l2dist(unsigned long x, unsigned long y);
  unsigned long long	l2dist(unsigned long long x, unsigned long long y);
  float			l2dist(float x, float y);
  double		l2dist(double x, double y);
  long double		l2dist(long double x, long double y);


  int			l2norm_sqr(int x);
  long			l2norm_sqr(long n);
  long long		l2norm_sqr(long long n);
  unsigned int		l2norm_sqr(unsigned int x);
  unsigned long		l2norm_sqr(unsigned long n);
  unsigned long long	l2norm_sqr(unsigned long long n);
  float			l2norm_sqr(float x);
  double		l2norm_sqr(double x);
  long double		l2norm_sqr(long double x);

  int			l2dist_sqr(int x, int y);
  long			l2dist_sqr(long x, long y);
  long long		l2dist_sqr(long long x, long long y);
  unsigned int		l2dist_sqr(unsigned int x, unsigned int y);
  unsigned long		l2dist_sqr(unsigned long x, unsigned long y);
  unsigned long long	l2dist_sqr(unsigned long long x, unsigned long long y);
  float			l2dist_sqr(float x, float y);
  double		l2dist_sqr(double x, double y);
  long double		l2dist_sqr(long double x, long double y);

  int			linfnorm(int x);
  long			linfnorm(long n);
  long long		linfnorm(long long n);
  unsigned int		linfnorm(unsigned int x);
  unsigned long		linfnorm(unsigned long n);
  unsigned long long	linfnorm(unsigned long long n);
  float			linfnorm(float x);
  double		linfnorm(double x);
  long double		linfnorm(long double x);

  template <unsigned p> int			lpnorm(int x);
  template <unsigned p> long			lpnorm(long x);
  template <unsigned p> long long		lpnorm(long long x);
  template <unsigned p> unsigned int		lpnorm(unsigned int x);
  template <unsigned p> unsigned long		lpnorm(unsigned long n);
  template <unsigned p> unsigned long long	lpnorm(unsigned long long n);
  template <unsigned p> float			lpnorm(float x);
  template <unsigned p> double			lpnorm(double x);
  template <unsigned p> long double		lpnorm(long double x);

  namespace functional
  {
    template <typename T = void> struct sum_t;
    template <typename T = void> struct prod_t;
    template <typename T = void> struct minimum_t;
    template <typename T = void> struct maximum_t;
    template <typename T = void> struct l0norm_t;
    template <typename T = void> struct l1norm_t;
    template <typename T = void> struct l2norm_t;
    template <typename T = void> struct l2norm_sqr_t;
    template <typename T = void> struct linfnorm_t;
    template <unsigned p, typename T = void> struct lpnorm_t;
  };


  /********************************/
  /*** Implementation            **/
  /********************************/

# define MLN_GEN_CODE(FUN, TYPE, OP)		\
  inline TYPE FUN(TYPE x) { return OP; }

# define MLN_GEN_BINARY_CODE(FUN, TYPE, OP)		\
  inline TYPE FUN(TYPE x, TYPE y) { return OP; }

# define MLN_GEN_CODE_2(TEMPLATE, FUN, TYPE, OP)	\
  TEMPLATE inline TYPE FUN(TYPE x) { return OP; }

# define MLN_FUNCTIONAL_GEN_CODE(FUN)                      \
  namespace functional {                                   \
  using mln::FUN;                                          \
  template <typename T>                                    \
  struct FUN##_t {                                         \
    typedef decltype(FUN(std::declval<T>())) result_type;  \
    auto operator() (const T& x) const -> decltype(FUN(x)) \
    { return FUN(x); }                                     \
  };                                                       \
  template <>                                              \
  struct FUN##_t<void> {                                   \
    template <typename T>                                  \
    auto operator() (const T& x) const -> decltype(FUN(x)) \
    { return FUN(x); }                                     \
  };                                                       \
  }

# define MLN_FUNCTIONAL_GEN_BINARY_CODE(FUN)               \
  namespace functional {                                   \
  using mln::FUN;                                          \
  template <typename T>                                    \
  struct FUN##_t {                                         \
    typedef decltype(FUN(std::declval<T>(), std::declval<T>())) result_type; \
    auto operator() (const T& x, const T& y) const -> decltype(FUN(x,y)) \
    { return FUN(x,y); }                                                \
  };                                                       \
  template <>                                              \
  struct FUN##_t<void> {                                   \
    template <typename T>                                  \
    auto operator() (const T& x, const T& y) const -> decltype(FUN(x,y)) \
    { return FUN(x,y); }                                                \
  };                                                       \
  }

# define MLN_FUNCTIONAL_GEN_CODE_2(TTYPE, TNAME, FUN)               \
  namespace functional {                                            \
    template <TTYPE TNAME, typename T>                              \
    struct FUN##_t {                                                \
      typedef decltype(FUN<TTYPE>(std::declval<T>())) result_type;  \
      auto operator() (const T& x) const -> decltype(FUN<TTYPE>(x)) \
      { return FUN<TTYPE>(x); }                                     \
    };                                                              \
    template <TTYPE TNAME>                                          \
    struct FUN##_t<TNAME, void>{                                    \
      template <typename T>                                         \
      auto operator() (const T& x) const -> decltype(FUN<TTYPE>(x)) \
      { return FUN<TTYPE>(x); }                                     \
    };                                                              \
  }

  MLN_FUNCTIONAL_GEN_CODE(sqrt);
  MLN_FUNCTIONAL_GEN_CODE(cbrt);
  MLN_FUNCTIONAL_GEN_CODE(pow);

  MLN_GEN_CODE(sqr, int, x*x);
  MLN_GEN_CODE(sqr, long, x*x);
  MLN_GEN_CODE(sqr, long long, x*x);
  MLN_GEN_CODE(sqr, unsigned int, x*x);
  MLN_GEN_CODE(sqr, unsigned long, x*x);
  MLN_GEN_CODE(sqr, unsigned long long, x*x);
  MLN_GEN_CODE(sqr, float, x*x);
  MLN_GEN_CODE(sqr, double, x*x);
  MLN_GEN_CODE(sqr, long double, x*x);
  MLN_FUNCTIONAL_GEN_CODE(sqr);

  MLN_GEN_CODE(abs, int, std::abs(x));
  MLN_GEN_CODE(abs, long, std::abs(x));
  MLN_GEN_CODE(abs, long long, std::abs(x));
  MLN_GEN_CODE(abs, unsigned int, x);
  MLN_GEN_CODE(abs, unsigned long, x);
  MLN_GEN_CODE(abs, unsigned long long, x);
  MLN_GEN_CODE(abs, float, std::abs(x));
  MLN_GEN_CODE(abs, double, std::abs(x));
  MLN_GEN_CODE(abs, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE(abs);

  MLN_GEN_CODE(l0norm, int, std::abs(x));
  MLN_GEN_CODE(l0norm, long, std::abs(x));
  MLN_GEN_CODE(l0norm, long long, std::abs(x));
  MLN_GEN_CODE(l0norm, unsigned int, x);
  MLN_GEN_CODE(l0norm, unsigned long, x);
  MLN_GEN_CODE(l0norm, unsigned long long, x);
  MLN_GEN_CODE(l0norm, float, std::abs(x));
  MLN_GEN_CODE(l0norm, double, std::abs(x));
  MLN_GEN_CODE(l0norm, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE(l0norm);

  MLN_GEN_CODE(l1norm, int, std::abs(x));
  MLN_GEN_CODE(l1norm, long, std::abs(x));
  MLN_GEN_CODE(l1norm, long long, std::abs(x));
  MLN_GEN_CODE(l1norm, unsigned int, x);
  MLN_GEN_CODE(l1norm, unsigned long, x);
  MLN_GEN_CODE(l1norm, unsigned long long, x);
  MLN_GEN_CODE(l1norm, float, std::abs(x));
  MLN_GEN_CODE(l1norm, double, std::abs(x));
  MLN_GEN_CODE(l1norm, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE(l1norm);

  MLN_GEN_CODE(l2norm, int, std::abs(x));
  MLN_GEN_CODE(l2norm, long, std::abs(x));
  MLN_GEN_CODE(l2norm, long long, std::abs(x));
  MLN_GEN_CODE(l2norm, unsigned int, x);
  MLN_GEN_CODE(l2norm, unsigned long, x);
  MLN_GEN_CODE(l2norm, unsigned long long, x);
  MLN_GEN_CODE(l2norm, float, std::abs(x));
  MLN_GEN_CODE(l2norm, double, std::abs(x));
  MLN_GEN_CODE(l2norm, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE(l2norm);

  MLN_GEN_BINARY_CODE(l2dist, int, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, long, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, long long, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, unsigned int, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, unsigned long, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, unsigned long long, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, float, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, double, abs(x - y));
  MLN_GEN_BINARY_CODE(l2dist, long double, abs(x - y));
  MLN_FUNCTIONAL_GEN_BINARY_CODE(l2dist);

  MLN_GEN_CODE(l2norm_sqr, int, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, long, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, long long, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, unsigned int, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, unsigned long, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, unsigned long long, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, float, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, double, sqr(x));
  MLN_GEN_CODE(l2norm_sqr, long double, sqr(x));
  MLN_FUNCTIONAL_GEN_CODE(l2norm_sqr);

  MLN_GEN_BINARY_CODE(l2dist_sqr, int, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, long, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, long long, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, unsigned int, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, unsigned long, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, unsigned long long, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, float, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, double, sqr(x - y));
  MLN_GEN_BINARY_CODE(l2dist_sqr, long double, sqr(x - y));
  MLN_FUNCTIONAL_GEN_BINARY_CODE(l2dist_sqr);

  MLN_GEN_CODE(linfnorm, int, std::abs(x));
  MLN_GEN_CODE(linfnorm, long, std::abs(x));
  MLN_GEN_CODE(linfnorm, long long, std::abs(x));
  MLN_GEN_CODE(linfnorm, unsigned int, x);
  MLN_GEN_CODE(linfnorm, unsigned long, x);
  MLN_GEN_CODE(linfnorm, unsigned long long, x);
  MLN_GEN_CODE(linfnorm, float, std::abs(x));
  MLN_GEN_CODE(linfnorm, double, std::abs(x));
  MLN_GEN_CODE(linfnorm, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE(linfnorm);

  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, int, std::abs(x));
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, long, std::abs(x));
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, long long, std::abs(x));
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, unsigned int, x);
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, unsigned long, x);
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, unsigned long long, x);
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, float, std::abs(x));
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, double, std::abs(x));
  MLN_GEN_CODE_2(template <unsigned p>, lpnorm, long double, std::abs(x));
  MLN_FUNCTIONAL_GEN_CODE_2(unsigned, p, lpnorm);

  MLN_GEN_CODE(maximum, int, x);
  MLN_GEN_CODE(maximum, long, x);
  MLN_GEN_CODE(maximum, long long, x);
  MLN_GEN_CODE(maximum, unsigned int, x);
  MLN_GEN_CODE(maximum, unsigned long, x);
  MLN_GEN_CODE(maximum, unsigned long long, x);
  MLN_GEN_CODE(maximum, float, x);
  MLN_GEN_CODE(maximum, double, x);
  MLN_GEN_CODE(maximum, long double, x);
  MLN_FUNCTIONAL_GEN_CODE(maximum);

  MLN_GEN_CODE(minimum, int, x);
  MLN_GEN_CODE(minimum, long, x);
  MLN_GEN_CODE(minimum, long long, x);
  MLN_GEN_CODE(minimum, unsigned int, x);
  MLN_GEN_CODE(minimum, unsigned long, x);
  MLN_GEN_CODE(minimum, unsigned long long, x);
  MLN_GEN_CODE(minimum, float, x);
  MLN_GEN_CODE(minimum, double, x);
  MLN_GEN_CODE(minimum, long double, x);
  MLN_FUNCTIONAL_GEN_CODE(minimum);

  MLN_GEN_CODE(sum, int, x);
  MLN_GEN_CODE(sum, long, x);
  MLN_GEN_CODE(sum, long long, x);
  MLN_GEN_CODE(sum, unsigned int, x);
  MLN_GEN_CODE(sum, unsigned long, x);
  MLN_GEN_CODE(sum, unsigned long long, x);
  MLN_GEN_CODE(sum, float, x);
  MLN_GEN_CODE(sum, double, x);
  MLN_GEN_CODE(sum, long double, x);
  MLN_FUNCTIONAL_GEN_CODE(sum);

  MLN_GEN_CODE(prod, int, x);
  MLN_GEN_CODE(prod, long, x);
  MLN_GEN_CODE(prod, long long, x);
  MLN_GEN_CODE(prod, unsigned int, x);
  MLN_GEN_CODE(prod, unsigned long, x);
  MLN_GEN_CODE(prod, unsigned long long, x);
  MLN_GEN_CODE(prod, float, x);
  MLN_GEN_CODE(prod, double, x);
  MLN_GEN_CODE(prod, long double, x);
  MLN_FUNCTIONAL_GEN_CODE(prod);

# undef MLN_GEN_CODE
# undef MLN_GEN_BINARY_CODE
# undef MLN_GEN_CODE_2
# undef MLN_FUNCTIONAL_GEN_CODE
# undef MLN_FUNCTIONAL_GEN_BINARY_CODE
# undef MLN_FUNCTIONAL_GEN_CODE_2
}

#endif // ! MLN_CORE_MATH_OPS_HPP
