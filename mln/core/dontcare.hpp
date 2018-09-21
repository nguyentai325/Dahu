#ifndef MLN_CORE_DONTCARE_HPP
# define MLN_CORE_DONTCARE_HPP

/// \file Dontcare type (a dummy type used in metaprogramming)

namespace mln
{

  /// \brief dontcare type.
  ///
  /// It can (implicitely) constucted from any type.
  struct dontcare_t
  {

    template <class... T>
    constexpr
    dontcare_t(T...)
    {
    };

    template <class... T>
    void
    operator() (T...) const
    {
    }

  };

  namespace
  {
    static constexpr dontcare_t dontcare = {};
  }

} // end of namespace mln

#endif //!MLN_CORE_DONTCARE_HPP
