#ifndef VALUE_TRAITS_HPP
# define VALUE_TRAITS_HPP

# include <functional>
# include <limits>

namespace mln
{
  // Forward declaration: special ordering.
  template <typename U, typename V = U>
  struct productorder_less;

  template <typename U, typename V = U>
  struct productorder_less_equal;

  template <typename U, typename V = U>
  struct lexicographicalorder_less;



  template <typename V, typename Ordering = productorder_less<V>,
            class Enable = void>
  struct value_traits
  {
    static_assert(!std::is_same<Enable, void>::value,
		  "You must specialize this trait for your types.");
  };

  // Default traits for std::greater
  template <typename V, class Enable>
  struct value_traits<V, std::greater<V>, Enable>
  {
    static constexpr unsigned ndim = value_traits<V, std::less<V>, Enable>::ndim;
    static constexpr unsigned quant = value_traits<V, std::less<V>, Enable>::quant;
    static constexpr V min() { return value_traits<V, std::less<V>, Enable>::max(); }
    static constexpr V max() { return value_traits<V, std::less<V>, Enable>::min(); }
    static constexpr V inf() { return value_traits<V, std::less<V>, Enable>::sup(); }
    static constexpr V sup() { return value_traits<V, std::less<V>, Enable>::inf(); }

  };


  // Specialization for primitive types
  template <typename V>
  struct value_traits<V, std::less<V>, typename std::enable_if<std::is_arithmetic<V>::value>::type>
  {
    static constexpr unsigned ndim = 1;
    static constexpr unsigned quant = sizeof(V) * 8;
    static constexpr V min() { return std::numeric_limits<V>::lowest(); }
    static constexpr V max() { return std::numeric_limits<V>::max(); }
    static constexpr V inf() { return min(); }
    static constexpr V sup() { return max(); }
  };

  template <typename V>
  struct value_traits<V, productorder_less<V>, typename std::enable_if<std::is_arithmetic<V>::value>::type>
    : value_traits<V, std::less<V> >
  {
  };

  template <typename V>
  struct value_traits<V, lexicographicalorder_less<V>, typename std::enable_if<std::is_arithmetic<V>::value>::type>
    : value_traits<V, std::less<V> >
  {
  };


}

#endif // ! VALUE_TRAITS_HPP
