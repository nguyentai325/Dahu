#ifndef MLN_CORE_ALWAYS_HPP
# define MLN_CORE_ALWAYS_HPP

namespace mln
{

  /// \brief Function object that always return the same value
  template <class T>
  struct always_t;

  /// \brief Function object that always return true
  struct yes_t;

  /// \brief Function object that always return false
  struct no_t;

  /*******************************/
  /***  Implementation         ***/
  /*******************************/

  template <class T>
  struct always_t
  {
    constexpr
    always_t(T x = T())
    : m_x(x)
    {
    }

    template <class... TArgs>
    constexpr
    T
    operator() (TArgs...) const
    {
      return m_x;
    }

  private:
    T m_x;
  };

  struct yes_t : always_t<bool>
  {
    constexpr yes_t()
      : always_t<bool>(true)
    {
    }
  };

  struct no_t : always_t<bool>
  {
    constexpr no_t()
      : always_t<bool>(false)
    {
    }
  };

  static constexpr no_t  no;
  static constexpr yes_t yes;

}

#endif // ! MLN_CORE_ALWAYS_HPP
