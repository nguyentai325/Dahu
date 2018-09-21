#ifndef MLN_CORE_IMAGE_DMORPH_EXTENDED_BY_VALUE_IMAGE_HPP
# define MLN_CORE_IMAGE_DMORPH_EXTENDED_BY_VALUE_IMAGE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image_base.hpp>

namespace mln
{

  template <typename I>
  struct extended_by_value_image;

  /******************************************/
  /****              Traits              ****/
  /******************************************/


  template <typename I>
  struct image_traits< extended_by_value_image<I> >
  {
    typedef typename image_traits<I>::category   category;
    typedef typename image_traits<I>::accessible accessible;
    typedef typename image_traits<I>::indexable  indexable;
    typedef std::false_type                      concrete;
    typedef extension::value_extension_tag       extension;
  };


  template <typename I>
  struct image_concrete< extended_by_value_image<I> >
  {
    typedef typename image_concrete<I>::type    type;
  };

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename I>
  struct extended_by_value_image : image_base< extended_by_value_image<I>,
                                               typename std::remove_reference<I>::type::point_type,
                                               typename std::remove_reference<I>::type::value_type>
  {
  private:
    typedef typename std::remove_reference<I>::type     image_t;

  public:
    typedef typename image_value<image_t>::type              value_type;
    typedef typename image_reference<image_t>::type          reference;
    typedef typename image_const_reference<image_t>::type    const_reference;
    typedef typename image_pixel<image_t>::type              pixel_type;
    typedef typename image_const_pixel<image_t>::type        const_pixel_type;

    typedef typename image_value_range<image_t>::type        value_range;
    typedef typename image_const_value_range<image_t>::type  const_value_range;
    typedef typename image_pixel_range<image_t>::type        pixel_range;
    typedef typename image_const_pixel_range<image_t>::type  const_pixel_range;





  };


} // end of namespace mln

#endif //!MLN_CORE_IMAGE_DMORPH_EXTENDED_BY_VALUE_IMAGE_HPP
