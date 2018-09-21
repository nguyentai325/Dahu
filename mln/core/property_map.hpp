#ifndef MLN_CORE_PROPERTY_MAP_HPP
# define MLN_CORE_PROPERTY_MAP_HPP

namespace mln
{

  template <class V, class>
  struct property_map
  {
    static_assert(not std::is_same<V, V>::value, "This class must be specialized.");
  };


  /// \brief An helper to convert a function \p f to a property map \p pm.
  /// It simply provides a functor with \p pm[x] that redirects to \p f(x)
  template <class F, class ArgType>
  struct functional_property_map
  {
    typedef typename std::result_of<F(ArgType)>::type                   reference;
    typedef typename std::result_of<const F(ArgType)>::type             const_reference;
    typedef typename std::decay<reference>::type                        value_type;
    typedef ArgType                                                     key_type;

    functional_property_map(const F& fun)
    : f (fun)
    {
    }

    template <class... T>
    typename std::result_of<const F(T...)>::type
    operator[] (T&&... args) const
    {
      return f(std::forward<T>(args)...);
    }

    template <class... T>
    typename std::result_of<F(T...)>::type
    operator[] (T&&... args)
    {
      return f(std::forward<T>(args)...);
    }

  private:
    F f;
  };

  template <class F>
  functional_property_map<F, typename F::argument_type>
  make_functional_property_map(const F& fun)
  {
    return functional_property_map<F, typename F::argument_type>(fun);
  }

  template <class ArgType, class F>
  functional_property_map<F, ArgType>
  make_functional_property_map(const F& fun)
  {
    return functional_property_map<F, ArgType>(fun);
  }


}

#endif // ! MLN_CORE_PROPERTY_MAP_HPP
