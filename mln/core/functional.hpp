#ifndef FUNCTIONAL_HPP
# define FUNCTIONAL_HPP

# include <functional>

namespace mln
{

  //
  // contrary to std::bind, it can be DefaultConstructible and Assignable
  template <typename BinaryFunction, typename V>
  struct binder1st
  {
  private:
    BinaryFunction fun_;
    V v_;

  public:
    binder1st() = default;
    binder1st(const BinaryFunction& f, const V& v)
      : fun_(f), v_ (v)
    {
    }

    template <typename T>
    typename std::result_of< BinaryFunction(V, T) >::type
    operator () (T&& x) const
    {
      return fun_(v_, std::forward<T>(x));
    }
  };

  template <typename BinaryFunction, typename V>
  struct binder2nd
  {
  private:
    BinaryFunction fun_;
    V v_;

  public:
    binder2nd() = default;

    binder2nd(const BinaryFunction& f, const V& v)
      : fun_(f), v_ (v)
    {
    }

    template <typename T>
    typename std::result_of< BinaryFunction(T, V) >::type
    operator () (T&& x) const
    {
      return fun_(std::forward<T>(x), v_);
    }
  };

  template <typename BinaryFunction, typename V>
  inline
  binder1st<BinaryFunction, V>
  bind1st(const BinaryFunction& f, const V& v)
  {
    return binder1st<BinaryFunction, V>(f, v);
  }

  template <typename BinaryFunction, typename V>
  inline
  binder2nd<BinaryFunction, V>
  bind2nd(const BinaryFunction& f, const V& v)
  {
    return binder2nd<BinaryFunction, V>(f, v);
  }

}

#endif // ! FUNCTIONAL_HPP
