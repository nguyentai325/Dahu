#ifndef MLN_CORE_FOREACH_HPP
# define MLN_CORE_FOREACH_HPP

# include <utility>
# include <mln/core/range/iter.hpp>

namespace mln
{

  namespace internal
  {

    template <typename T>
    struct false_var_t
    {
      false_var_t(const T& v) : x_(v) {}

      constexpr operator bool() const { return false; }
      T& get() { return x_;}
      false_var_t& set(const T& v) { x_ = v; return *this; }

      T x_;
    };

    template <typename ColExpr>
    struct should_copy_col_;

    // The expression is a lvalue, no need to copy
    template <typename ColExpr>
    struct should_copy_col_<ColExpr&>
    {
      typedef ColExpr& type;

      static ColExpr&
      copy(ColExpr& col) { return col; }
    };

    // The expression is a rvalue, need to copy
    template <typename ColExpr>
    struct should_copy_col_
    {
      typedef ColExpr type;

      static ColExpr
        copy(ColExpr& col) { return std::move(col); }
    };

    template <typename ColExpr>
    typename should_copy_col_<ColExpr>::type
    should_copy_col(ColExpr&& x)
    {
      return should_copy_col_<ColExpr>::copy(x);
    }

  }

}

# define MLN_DECL_VAR(ID, VALUE)				\
  if (mln::internal::false_var_t< decltype(VALUE) > ID = VALUE) {} else


# define __mln_should_copy_col__(COL, ID)          \
  decltype(mln::internal::should_copy_col(COL)) ID = mln::internal::should_copy_col(COL)

# define __mln_should_copy_col_local__(COL, ID)				\
  if (mln::internal::false_var_t< decltype(mln::internal::should_copy_col(COL)) > \
      ID = mln::internal::should_copy_col(COL) ) {} else


/******************************************/
/****         mln_foreach macro         ****/
/******************************************/

# define __mln_do_local__(EXPR)			\
  if ((EXPR), false) {} else


# define mln_foreach(p, COL)						\
  __mln_should_copy_col_local__(COL, _mln_range_)			\
  MLN_DECL_VAR(_mln_it_, mln::rng::iter(_mln_range_.get()))		\
  MLN_DECL_VAR(_mln_continue_, true)					\
  for (_mln_it_.get().init();						\
       _mln_continue_.get() and !_mln_it_.get().finished();		\
       _mln_continue_.get() ? _mln_it_.get().next() : (void) 0)		\
    if (_mln_continue_.set(false)) {} else				\
      for (p = *(_mln_it_.get()); !_mln_continue_.get(); _mln_continue_.set(true)) \


# define mln_reverse_foreach(p, COL)						\
  __mln_should_copy_col_local__(COL, _mln_range_)			\
  MLN_DECL_VAR(_mln_it_, mln::rng::riter(_mln_range_.get()))		\
  MLN_DECL_VAR(_mln_continue_, true)					\
  for (_mln_it_.get().init();						\
       _mln_continue_.get() and !_mln_it_.get().finished();		\
       _mln_continue_.get() ? _mln_it_.get().next() : (void) 0)		\
    if (_mln_continue_.set(false)) {} else				\
      for (p = *(_mln_it_.get()); !_mln_continue_.get(); _mln_continue_.set(true))


#endif // ! MLN_CORE_FOREACH_HPP
