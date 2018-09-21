#ifndef MLN_CORE_OBJECT_WRAPPERS_HPP
# define MLN_CORE_OBJECT_WRAPPERS_HPP

// Uniform object holders to mimics reference wrapper.

namespace mln
{

  template <class T>
  struct object_wrapper;


  /****************************/
  /*** Implementation       ***/
  /****************************/

  template <class T>
  struct object_wrapper
  {
    typedef T type;

    object_wrapper() = default;

    object_wrapper(const T& x)
    : m_x (x)
    {
    }

    object_wrapper(T&& x)
    : m_x (std::move(x))
    {
    }

    object_wrapper&
    operator= (const T& x)
    {
      m_x = x;
      return *this;
    }

    object_wrapper&
    operator= (T&& x)
    {
      m_x = std::move(x);
      return *this;
    }

    operator T& ()
    {
      return m_x;
    }

    operator const T& () const
    {
      return m_x;
    }

    T& get()
    {
      return m_x;
    }

    const T& get() const
    {
      return m_x;
    }

    T m_x;
  };

}

#endif // ! MLN_CORE_OBJECT_WRAPPERS_HPP
