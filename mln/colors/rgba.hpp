#ifndef MLN_COLORS_RGBA_HPP
# define MLN_COLORS_RGBA_HPP

# include <mln/core/value/int.hpp>
# include <mln/core/colors.hpp>
# include <mln/core/ops.hpp>
# include <mln/core/vec_base.hpp>
# include <mln/core/image/morphers/transformed_image.hpp>

namespace mln
{

  namespace colors
  {
    struct rgba_tag {};
  }

  namespace internal
  {
    template <>
    struct vec_base_traits<mln::colors::rgba_tag>
      {
        static const bool is_additive = true;
        static const bool is_additive_ext = true;
        static const bool is_multiplicative = false;
        static const bool is_multiplicative_ext = true;
        static const bool is_less_than_comparable = false;
        static const bool is_equality_comparable = true;
    };

  }

  namespace colors
  {
    template <typename T>
    using rgba = internal::vec_base<T, 4, rgba_tag>;

    typedef rgba<uint8> rgba8;
    typedef rgba<uint16> rgba16;

  }

}


MLN_DECLARE_IMAGE_LVALUE_OPERATOR_OVERLOAD
(red, (mln::colors::rgba8), (mln::getter<0>))

MLN_DECLARE_IMAGE_LVALUE_OPERATOR_OVERLOAD
(green, (mln::colors::rgba8), (mln::getter<1>))

MLN_DECLARE_IMAGE_LVALUE_OPERATOR_OVERLOAD
(blue, (mln::colors::rgba8), (mln::getter<2>))


#endif // ! MLN_COLORS_RGBA_HPP
