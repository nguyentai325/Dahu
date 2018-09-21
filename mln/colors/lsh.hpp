#ifndef MLN_COLORS_LSH_HPP
# define MLN_COLORS_LSH_HPP

# include <mln/core/vec_base.hpp>
# include <mln/core/colors.hpp>

// FIXME: optimize this out (slow because of floats and saturations)

namespace mln
{

  struct lsh_tag {};

  template <typename T>
  using lsh = internal::vec_base<T, 3, lsh_tag>;

  typedef lsh<uint8> lsh8;

  namespace internal
  {
    template <>
    struct vec_base_traits<lsh_tag>
    {
      static const bool is_additive = true;
      static const bool is_additive_ext = true;
      static const bool is_multiplicative = false;
      static const bool is_multiplicative_ext = true;
      static const bool is_less_than_comparable = true;
      static const bool is_equality_comparable = true;
    };

  }

  template <typename T>
  lsh<T> rgb2lsh(const rgb<T>& v);

  template <typename T>
  rgb<T> lsh2rgb(const lsh<T>& v);


  /*********************/
  /*** Implementation **/
  /*********************/

  namespace internal
  {
    template <typename T>
    inline
    std::tuple<T,T,T,int>
    get_MinMedMax(T r, T g, T b)
    {
      if (r < g) {
	if (g < b)	return std::make_tuple(r, g, b, 3); // r < g < b
	else if (b < r)	return std::make_tuple(b, r, g, 1); // b < r < g
	else		return std::make_tuple(r, b, g, 2); // r < b <= g ou r <= b < g
      } else { // g <= r
	if (b < g)	return std::make_tuple(b, g, r, 0); // b < g <= r
	else if (r < b) return std::make_tuple(g, r, b, 4); // g <= r < b
	else		return std::make_tuple(g, b, r, 5); // g <= b <= r
      }
    }
  }

  template <typename T>
  inline
  lsh<T>
  rgb2lsh(const rgb<T>& v)
  {
    static constexpr float k = (1 << value_traits<T>::quant) / 6;
    static constexpr int   Hnv = (int)k * 6;

    T min, med, max;
    int lambda;
    std::tie(min, med, max, lambda) = internal::get_MinMedMax(v[0],v[1],v[2]);
    float L = (min + max + med) / 3.0;
    float S = 1.5f * (L > med ? (max-L) : (L-min));
    float H = k * (lambda + 0.5 - ((lambda % 2 == 0) ? 1 : -1) * (max + min - 2*med) / (2*S));

    lsh<T> out;
    out[0] = L;
    out[1] = std::floor(S);
    out[2] = (int)std::floor(H) % Hnv;
    return out;
  }


  template <typename T>
  inline
  rgb<T>
  lsh2rgb(const lsh<T>& v)
  {
    static constexpr float k = (1 << value_traits<T>::quant) / 6;

    rgb<T> out;

    T l = v[0], s = v[1], h = v[2];

    int lambda = h / (int) k;
    int m = lambda % 2 == 0 ? 1 : -1;
    float phi = h / k - lambda;
    float s1_3 = s * (1.0/3), s2_3 = s * (2.0/3);
    float c = m * (phi - 0.5);

    float maxi_ = (c <= 0) ? (l + s2_3) : (l + s1_3 + s2_3 * (((lambda+1)%2) - m * phi));
    float mini_ = (c >= 0) ? (l - s2_3) : (l - s1_3 - s2_3 * (lambda%2 + m * phi));
    float med_ = l - m * s1_3 + m * s2_3 * phi;

    T max = std::max(0, std::min(255, (int) (maxi_ + 0.5) ));
    T med = std::max(0, std::min(255, (int) (med_ + 0.5) ));
    T min = std::max(0, std::min(255, (int) (mini_ + 0.5) ));

    switch (lambda) {
      case 0: return rgb<T> {max,med,min};
      case 1: return rgb<T> {med,max,min};
      case 2: return rgb<T> {min,max,med};
      case 3: return rgb<T> {min,med,max};
      case 4: return rgb<T> {med,min,max};
      case 5: return rgb<T> {max,min,med};
      default: assert(false); return rgb<T> ();
    }
  }

}

#endif // ! MLN_COLORS_LSH_HPP
