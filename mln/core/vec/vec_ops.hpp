#ifndef MLN_CORE_VEC_VEC_OPS_HPP
# define MLN_CORE_VEC_VEC_OPS_HPP

namespace mln
{

  struct dyn_getter
  {
    dyn_getter(int k) : m_k(k) {}

    template <class T, unsigned dim, class Tag>
    const T&
    operator () (const internal::vec_base<T, dim, Tag>& x) const
    {
      return x[m_k];
    }

    template <class T, unsigned dim, class Tag>
    T&
    operator () (internal::vec_base<T, dim, Tag>& x) const
    {
      return x[m_k];
    }

    template <class T, unsigned dim, class Tag>
    T&&
    operator () (internal::vec_base<T, dim, Tag>&& x) const
    {
      return std::forward<T&&>(x[m_k]);
    }

    int m_k;
  };

}

#endif // ! MLN_CORE_VEC_VEC_OPS_HPP
