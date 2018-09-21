#ifndef MLN_CORE_ITERATOR_IMAGE_ACCESS_ITERATOR_HPP
# define MLN_CORE_ITERATOR_IMAGE_ACCESS_ITERATOR_HPP


/// \file
/// \ingroup iterators
/// \ingroup image_iterators
///
/// \brief Provide value iterator and pixel iterator
/// over an image from a custom domain iterator.
///
/// From a custom site iterator \p p, and an image \p ima
/// these classes provides:
/// * value iterator that maps to \p ima(p)
/// * a pixel iterator that maps to \p (p,ima(p))
///
/// TODO: add interopabilility

# include <mln/core/iterator/iterator_base.hpp>
# include <mln/core/image_traits.hpp>

namespace mln
{

  template <typename Image, typename SiteIterator>
  struct image_access_value_iterator;

  template <typename Image, typename SiteIterator, typename MorpherImage>
  struct image_access_pixel_iterator;


  /********************/
  /* Implementations  */
  /********************/

  template <typename Image, typename SiteIterator>
  struct image_access_value_iterator :
    iterator_base< image_access_value_iterator<Image, SiteIterator>,
		   typename image_value<Image>::type,
		   typename image_reference<Image>::type
		   >
  {
  public:
    typedef typename image_reference<Image>::type reference;

  public:
    image_access_value_iterator() = default;
    image_access_value_iterator(const image_access_value_iterator&) = default;

    image_access_value_iterator(Image& ima, SiteIterator p)
      : m_ima(&ima), m_p (p)
    {
    }

    template <class Image2, class SiteIterator2>
    image_access_value_iterator(const image_access_value_iterator<Image2, SiteIterator2>& other,
				typename std::enable_if< std::is_convertible<Image2*, Image*>::value and
							 std::is_convertible<SiteIterator2, SiteIterator>::value >::type* = NULL)
      : m_ima (other.m_ima), m_p (other.m_p)
    {
    }

    void init() { m_p.init(); }
    void next() { m_p.next(); }
    bool finished() const { return m_p.finished(); }
    reference dereference() const { return (*m_ima)(*m_p); }

  private:
    template <class, class> friend struct image_access_value_iterator;

    Image*	 m_ima;
    SiteIterator m_p;
  };


  template <typename Image, typename SiteIterator, typename MorpherImage>
  struct image_access_pixel
  {
    typedef typename image_value<Image>::type		value_type;
    typedef typename image_reference<Image>::type	reference;
    typedef typename Image::point_type			point_type;
    typedef typename Image::site_type			site_type;
    typedef MorpherImage				image_type;

    image_access_pixel() = default;
    image_access_pixel(const image_access_pixel&) = default;

    image_access_pixel(Image& ima, const SiteIterator& p, MorpherImage& morph)
      : m_morpher(&morph), m_ima(&ima), m_p(p)
    {
    }

    template <typename Image2, typename SiteIterator2, typename MorpherImage2>
    image_access_pixel(const image_access_pixel<Image2, SiteIterator2, MorpherImage2>& other,
		       typename std::enable_if<std::is_convertible<Image2*, Image*>::value and
		       std::is_convertible<SiteIterator2, SiteIterator>::value and
		       std::is_convertible<MorpherImage2*, MorpherImage*>::value>::type* = NULL)
      : m_morpher (other.m_morpher), m_ima (other.m_ima), m_p (other.m_p)
    {
    }


    reference   val() const { return (*m_ima)(*m_p); }
    site_type   point() const { return *m_p; }
    site_type   site() const { return *m_p; }
    image_type&  image() const { return *m_morpher; }

    template <class, class, class> friend struct image_access_pixel_iterator;
    template <class, class, class> friend struct image_access_pixel;

  private:
    MorpherImage* m_morpher;
    Image*	  m_ima;
    SiteIterator  m_p;
  };

  template <typename Image, typename SiteIterator, typename MorpherImage>
  struct image_access_pixel_iterator :
    iterator_base< image_access_pixel_iterator<Image, SiteIterator, MorpherImage>,
		   const image_access_pixel<Image, SiteIterator, MorpherImage> >
  {
  private:
    typedef image_access_pixel<Image, SiteIterator, MorpherImage> pixel_type;

  public:

    image_access_pixel_iterator() = default;
    image_access_pixel_iterator(Image& ima, const SiteIterator& p, MorpherImage& morpher)
      : m_pix(ima, p, morpher)
    {
    }

    template <class Image2, class SiteIterator2, class MorpherImage2>
    image_access_pixel_iterator(const image_access_pixel_iterator<Image2, SiteIterator2, MorpherImage2>& other)
      : m_pix (other.m_pix)
    {
    }

    void init() { m_pix.m_p.init(); }
    void next() { m_pix.m_p.next(); }
    bool finished() const { return m_pix.m_p.finished(); }
    const pixel_type& dereference() const { return m_pix; }

  private:
    template <class, class, class> friend struct image_access_pixel_iterator;

    pixel_type m_pix;
  };


}

#endif // ! MLN_CORE_ITERATOR_IMAGE_ACCESS_ITERATOR_HPP
