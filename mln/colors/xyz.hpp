#ifndef MLN_COLORS_XYZ_HPP
# define MLN_COLORS_XYZ_HPP

# include <mln/core/vec_base.hpp>
# include <mln/core/colors.hpp>

// FIXME: optimize this out (slow because of floats and saturations)

namespace mln
{

  struct xyz_tag {};

  template <typename T>
  using xyz = internal::vec_base<T, 3, xyz_tag>;

  typedef xyz<uint8> xyz8;

  namespace internal
  {
    template <>
    struct vec_base_traits<xyz_tag>
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
  xyz<float> rgb2xyz(const rgb<T>& v);

  template <typename T>
  rgb<float> xyz2rgb(const xyz<T>& v);


  /*********************/
  /*** Implementation **/
  /*********************/

  template <typename T>
  inline
  xyz<float>
  rgb2xyz(const rgb<T>& v)
  {
    float x = (0.4125 * v[0] + 0.3576 * v[1] + 0.1805 * v[2]);
    float y = (0.2127 * v[0] + 0.7152 * v[1] + 0.0722 * v[2]);
    float z = (0.0193 * v[0] + 0.1192 * v[1] + 0.9505 * v[2]);

    return {x,y,z};
  }

  template <typename T>
  inline
  rgb<float>
  xyz2rgb(const xyz<T>& v)
  {
    float r = ( 3.2406 * v[0] - 1.5372 * v[1] - 0.4986 * v[2]);
    float g = (-0.9689 * v[0] + 1.8758 * v[1] + 0.0415 * v[2]);
    float b = ( 0.0557 * v[0] - 0.2040 * v[1] + 1.0570 * v[2]);

    return {r,g,b};
  }

}

#endif // ! MLN_COLORS_XYZ_HPP
