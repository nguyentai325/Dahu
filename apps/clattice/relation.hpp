#ifndef RELATION_HPP
# define RELATION_HPP

# include <mln/core/vec.hpp>


/// \file Partial order definition on vectorial data type
/// using the product order.
/// u <= v      <=>      u[i] <= v[i] forall i
/// u < v       <=>      u <= v and u != v (by reflexive reduction)
///
/// Note that this is different from:
/// u < v       <=>      u[i] < v[i] forall i
/// u <= v      <=>      u < v or u = v (by reflexive closure)
/// With the second relation the sup/inf are not part of the set.

template <class Vec>
struct rng_porder_less
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = "<";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_isless(rng.upper, v);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return mln::vecprod_isless(u, v);
  }

  Vec inf(const Vec& u, const Vec& v) const
  {
    return mln::inf(u,v, mln::productorder_less<Vec> ());
  }
};

template <class Vec>
struct rng_porder_greater
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = ">";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_isgreater(rng.lower, v);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return mln::vecprod_isgreater(u, v);
  }

  Vec inf(const Vec& u, const Vec& v) const
  {
    return mln::sup(u,v, mln::productorder_less<Vec> ());
  }
};

template <class Vec>
struct rng_porder_less_equal
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = "<=";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_islessequal(rng.upper, v);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return mln::vecprod_islessequal(u, v);
  }

  Vec inf(const Vec& u, const Vec& v) const
  {
    return mln::inf(u,v, mln::productorder_less<Vec> ());
  }
};

template <class Vec>
struct rng_porder_greater_equal
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = ">=";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_isgreaterequal(rng.lower, v);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return mln::vecprod_isgreaterequal(u, v);
  }

  Vec inf(const Vec& u, const Vec& v) const
  {
    return mln::sup(u,v, mln::productorder_less<Vec> ());
  }
};

template <class Vec>
struct rng_porder_equal_to
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = "==";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_islessequal(rng.lower, v) and
      mln::vecprod_islessequal(v, rng.upper);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return u == v;
  }

  Vec inf(const Vec& u, const Vec&) const
  {
    return u;
  }

};

template <class Vec>
struct rng_porder_not_equal_to
{
  typedef mln::morpho::tos::irange<Vec> range_t;
  static constexpr const char* str = "!=";

  bool
  operator() (const range_t& rng, const Vec& v) const
  {
    return mln::vecprod_isless(v, rng.lower) or
      mln::vecprod_isgreater(v, rng.upper);
  }

  bool
  operator() (const Vec& u, const Vec& v) const
  {
    return u != v;
  }

  Vec inf(const Vec& u, const Vec& v) const
  {
    return u; // FIXME: false because != non transitive
  }

};

#endif // ! RELATION_HPP
