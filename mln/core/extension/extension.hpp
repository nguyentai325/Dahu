#ifndef MLN_CORE_EXTENSION_EXTENSION_HPP
# define MLN_CORE_EXTENSION_EXTENSION_HPP

# include <mln/core/extension/extension_traits.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/core/image/morphers/extended_by_value_image.hpp>
# include <mln/core/internal/get_border_from_nbh.hpp>

namespace mln
{

  namespace extension
  {

    /******************************************/
    /****          Free functions          ****/
    /******************************************/


    /// \brief Remove the extension of an image.
    ///
    /// \p remove_extension recursively removes the extensions of an image until
    /// getting an image without extension or an image whose extension cannot be
    /// removed. This function has to be overloaded by the morphers extending an
    /// image. This default implementation returns the input image as such.
    ///
    /// \param ima input image
    ///
    /// \return The image without extension
    ///
    template <typename I>
    const I&
    remove_extension(const Image<I>& ima);

    template <typename I>
    I&
    remove_extension(Image<I>& ima);

    template <typename I>
    I&&
    remove_extension(Image<I>&& ima);

    /// \brief Check if an image extension is wide enough to support
    /// a given neighborhood/se/window.
    ///
    /// An image does not need "adjutment" w.r.t to a neighborhood if
    /// these conditions are met:
    /// * \p ima has an extension
    /// * \p the neighborhood is constant (either static or dynamic)
    /// * \p the extension is wide enough
    template <class I, class N>
    bool need_adjust(const Image<I>& ima, const Neighborhood<N>& nbh);


    /// \brief Add an infinite value extension to the image. In the resulting
    /// image, every access outside the image domain yields in the extension value.
    ///
    /// \param ima input image
    /// \param value_extension_tag
    ///
    /// \return an image extended by value
    ///
    template <typename I>
    extended_by_value_image<const I&>
    add_value_extension(const Image<I>& ima, const mln_value(I)& v);

    template <typename I>
    extended_by_value_image<I&>
    add_value_extension(Image<I>& ima, const mln_value(I)& v);

    template <typename I>
    extended_by_value_image<I>
    add_value_extension(Image<I>&& ima, const mln_value(I)& v);


    /******************************************/
    /****          Implementation          ****/
    /******************************************/


    template <typename I>
    const I&
    remove_extension(const Image<I>& ima)
    {
      return exact(ima);
    }

    template <typename I>
    I&
    remove_extension(Image<I>& ima)
    {
      return exact(ima);
    }

    template <typename I>
    I&&
    remove_extension(Image<I>&& ima)
    {
      return move_exact(ima);
    }

    template <typename I>
    extended_by_value_image<const I&>
    add_value_extension(const Image<I>& ima, const mln_value(I)& v)
    {
      return extended_by_value_image<const I&>(exact(ima), v);
    }

    template <typename I>
    extended_by_value_image<I&>
    add_value_extension(Image<I>& ima, const mln_value(I)& v)
    {
      return extended_by_value_image<I&>(exact(ima), v);
    }

    template <typename I>
    extended_by_value_image<I>
    add_value_extension(Image<I>&& ima, const mln_value(I)& v)
    {
      return extended_by_value_image<I>(move_exact(ima), v);
    }

    namespace impl
    {

      template <class I, class N>
      bool need_adjust(const I&, const N&, extension_tag, adaptative_neighborhood_tag)
      {
        return true;
      }

      template <class I, class N>
      bool need_adjust(const I&, const N&, value_extension_tag, adaptative_neighborhood_tag)
      {
        return false;
      }

      template <class I, class N>
      bool need_adjust(const I& ima, const N& nbh, border_extension_tag, dynamic_neighborhood_tag)
      {
        return (int) internal::get_border_from_nbh(nbh) > (int) ima.border();
      }

    }


    template <class I, class N>
    bool need_adjust(const Image<I>& ima, const Neighborhood<N>& nbh)
    {
      return extension::impl::need_adjust(exact(ima), exact(nbh),
                                          typename image_traits<I>::extension (),
                                          typename N::category ());
    }

  } // end of namespace mln::extension
} // end of namespace mln

#endif //!MLN_CORE_EXTENSION_EXTENSION_HPP
