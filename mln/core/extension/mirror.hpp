#ifndef MLN_CORE_EXTENSION_MIRROR_HPP
# define MLN_CORE_EXTENSION_MIRROR_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/extension/extension_traits.hpp>
# include <mln/core/trace.hpp>
namespace mln
{
  namespace extension
  {

    template <typename I>
    const I& mirror(const Image<I>& ima);

    template <typename I>
    I&& mirror(Image<I>&& ima);

/******************************************/
/****          Implementation          ****/
/******************************************/


    template <typename I>
    const I&
    mirror(const Image<I>& ima)
    {
      static_assert(image_has_extension<I>::value, "Image must have an extension.");
      static_assert(extension_traits<typename image_extension_type<I>::type>::support_fill::value,
                    "Image extension must support filling.");
      mln_entering("mln::extension::mirror");
      exact(ima).extension().mirror();
      mln_exiting();
      return exact(ima);
    }

    template <typename I>
    I&&
    mirror(Image<I>&& ima)
    {
      static_assert(image_has_extension<I>::value, "Image must have an extension.");
      static_assert(extension_traits<typename image_extension_type<I>::type>::support_fill::value,
                    "Image extension must support mirror.");

      mln_entering("mln::extension::mirror");
      exact(ima).extension().mirror();
      mln_exiting();
      return move_exact(ima);
    }



  } // end of namespace mln::extension
} // end of namespace mln


#endif // ! MLN_CORE_EXTENSION_FILL_HPP
