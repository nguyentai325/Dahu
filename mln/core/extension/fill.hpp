#ifndef MLN_CORE_EXTENSION_FILL_HPP
# define MLN_CORE_EXTENSION_FILL_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/extension/extension_traits.hpp>

namespace mln
{
  namespace extension
  {

    template <typename I>
    const I& fill(const Image<I>& ima, mln_value(I) v);

    template <typename I>
    I&& fill(Image<I>&& ima, mln_value(I) v);

/******************************************/
/****          Implementation          ****/
/******************************************/


    template <typename I>
    const I&
    fill(const Image<I>& ima, mln_value(I) v)
    {
      static_assert(image_has_extension<I>::value, "Image must have an extension.");
      static_assert(extension_traits<typename image_extension_type<I>::type>::support_fill::value,
                    "Image extension must support filling.");

      exact(ima).extension().fill(v);
      return exact(ima);
    }

    template <typename I>
    I&&
    fill(Image<I>&& ima, mln_value(I) v)
    {
      static_assert(image_has_extension<I>::value, "Image must have an extension.");
      static_assert(extension_traits<typename image_extension_type<I>::type>::support_fill::value,
                    "Image extension must support filling.");

      exact(ima).extension().fill(v);
      return move_exact(ima);
    }



  } // end of namespace mln::extension
} // end of namespace mln


#endif // ! MLN_CORE_EXTENSION_FILL_HPP
