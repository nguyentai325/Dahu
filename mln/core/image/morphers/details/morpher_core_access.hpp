#ifndef MLN_CORE_IMAGE_MORPHERS_DETAILS_MORPHER_CORE_ACCESS_HPP
# define MLN_CORE_IMAGE_MORPHERS_DETAILS_MORPHER_CORE_ACCESS_HPP

namespace mln
{

    struct morpher_core_access
  {
    template <typename Morpher>
    static
    typename Morpher::image_t&
    get_ima_(Morpher* morpher)
    {
      return reinterpret_cast<typename Morpher::derived_t*>(morpher)->m_ima;
    }


    template <typename Morpher>
    static
    const typename Morpher::image_t&
    get_ima_(const Morpher* morpher)
    {
      return reinterpret_cast<const typename Morpher::derived_t*>(morpher)->m_ima;
    }

    template <typename Morpher>
    static
    typename Morpher::image_t&
    get_ima(Morpher* morpher)
    {
      return reinterpret_cast<typename Morpher::derived_t*>(morpher)->get_morphed();
    }


    template <typename Morpher>
    static
    const typename Morpher::image_t&
    get_ima(const Morpher* morpher)
    {
      return reinterpret_cast<const typename Morpher::derived_t*>(morpher)->get_morphed();
    }


    template <typename PixMorpher>
    static
    typename PixMorpher::pixel_t&
    get_pix_(PixMorpher* morpher)
    {
      return reinterpret_cast<typename PixMorpher::derived_t*>(morpher)->m_pix;
    }


    template <typename PixMorpher>
    static
    const typename PixMorpher::pixel_t&
    get_pix_(const PixMorpher* morpher)
    {
      return reinterpret_cast<const typename PixMorpher::derived_t*>(morpher)->m_pix;
    }

    template <typename PixMorpher>
    static
    typename PixMorpher::pixel_t&
    get_pix(PixMorpher* morpher)
    {
      return reinterpret_cast<typename PixMorpher::derived_t*>(morpher)->get_morphed();
    }


    template <typename PixMorpher>
    static
    const typename PixMorpher::pixel_t&
    get_pix(const PixMorpher* morpher)
    {
      return reinterpret_cast<const typename PixMorpher::derived_t*>(morpher)->get_morphed();
    }

  };


} // end of namespace mln

#endif //!MLN_CORE_IMAGE_MORPHERS_DETAILS_MORPHER_CORE_ACCESS_HPP
