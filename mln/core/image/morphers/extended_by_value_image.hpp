#ifndef EXTENDED_BY_VALUE_IMAGE_HPP
# define EXTENDED_BY_VALUE_IMAGE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/morphers/morpher_base.hpp>
# include <mln/core/range/iterator_range.hpp>

namespace mln
{


  /// \brief An image with an extension defined by a value
  template <typename I>
  struct extended_by_value_image;

  /// \{

  /// \brief Add an extension by value to an image
  template <typename I>
  extended_by_value_image<const I&>
  extend_by_value(const Image<I>& ima, const mln_value(I)& v);

  template <typename I>
  extended_by_value_image<I&>
  extend_by_value(Image<I>& ima, const mln_value(I)& v);

  template <typename I>
  extended_by_value_image<I>
  extend_by_value(Image<I>&& ima, const mln_value(I)& v);
  /// \}



  template <typename V>
  struct by_value_extension;


  /*****************************/
  /***  Traits               ***/
  /*****************************/

  template <typename I>
  struct image_traits< extended_by_value_image<I> >
  {
  private:
    typedef typename std::remove_reference<I>::type image_t;

  public:
    // Actually if I is concrete we could set this image concrete as well
    // but since adding this extension prevents indexation we better use
    // the concrete image of I rather than this morpher, as long as concretization
    // does not impose any property about this extension.

    typedef std::false_type				concrete;
    typedef typename std::common_type<typename image_traits<image_t>::category,
                                      bidirectional_image_tag>::type category;
    typedef typename image_traits<image_t>::accessible	accessible;

    // we can check if an index belongs to a domain using ima.domain().has(ima.point_at_index(i))
    // but this is not efficient. We better use points directly.
    typedef std::false_type				indexable;
    typedef mln::extension::value_extension_tag		extension;
  };


  template <typename I>
  struct image_concrete< extended_by_value_image<I> >
  {
    typedef mln_concrete(I) type;
  };

  namespace internal
  {
    template <typename I>
    struct image_init_from< extended_by_value_image<I> >
    {
      typedef typename image_init_from<typename std::decay<I>::type>::type type;
    };
  }


  /*****************************/
  /**** Implementation       ***/
  /*****************************/


  /******************************************/
  /****              Pixel              *****/
  /******************************************/

  template <typename Morpher, typename Pix>
  struct extended_by_value_image_pixel : Pixel< extended_by_value_image_pixel<Morpher, Pix> >
  {
    template <typename, typename>
    friend struct extended_by_value_image_pixel;

  public:
    typedef typename Pix::value_type value_type;
    typedef typename Pix::reference  reference;
    typedef typename Pix::point_type point_type;
    typedef typename Pix::site_type  site_type;
    typedef Morpher                  image_type;

    extended_by_value_image_pixel() = default;
    extended_by_value_image_pixel(const extended_by_value_image_pixel&) = default;

    extended_by_value_image_pixel(Morpher& ima, const Pix& pix)
      : m_ima(&ima), m_pix(pix)
    {
    }

    // Conversion non-const -> const pixel
    template <typename Morpher2, typename Pix2>
    extended_by_value_image_pixel(const extended_by_value_image_pixel<Morpher2, Pix2>& other,
                                  typename std::enable_if< std::is_convertible<Morpher2*, Morpher*>::value >::type* = NULL)
      : m_ima(other.m_ima), m_pix(other.m_pix)
    {
    }

    reference   val() const   { return m_pix.val(); }
    point_type  point() const { return m_pix.point(); }
    site_type   site() const  { return m_pix.site(); }
    image_type& image() const { return *m_ima; }

  private:
    Morpher*   m_ima;
    Pix        m_pix;
  };


  /******************************************/
  /****          Pixel iterators         ****/
  /******************************************/

  template <typename Morpher, typename PixelIterator>
  struct extended_by_value_image_pixel_iterator :
    iterator_base< extended_by_value_image_pixel_iterator<Morpher, PixelIterator>,
                   extended_by_value_image_pixel<Morpher, typename PixelIterator::value_type>,
                   extended_by_value_image_pixel<Morpher, typename PixelIterator::value_type>
                   >
  {
  private:
    typedef extended_by_value_image_pixel<Morpher, typename PixelIterator::value_type> pixel_type;
    template <typename, typename>
    friend struct extended_by_value_image_pixel_iterator;

  public:
    typedef typename PixelIterator::has_NL has_NL;

    extended_by_value_image_pixel_iterator() = default;
    extended_by_value_image_pixel_iterator(const extended_by_value_image_pixel_iterator&) = default;

    extended_by_value_image_pixel_iterator(Morpher& ima, const PixelIterator& it)
      : m_ima(&ima), m_it (it)
    {
    }

    template <typename Morpher2, typename PixelIterator2>
    extended_by_value_image_pixel_iterator(const extended_by_value_image_pixel_iterator<Morpher2, PixelIterator2>& other,
                                           typename std::enable_if< std::is_convertible<Morpher2*, Morpher*>::value>::type* = NULL)
      : m_ima(other.m_ima), m_it (other.m_it)
    {
    }


    void init()           { m_it.init(); }
    void next()           { m_it.next(); }
    bool finished() const { return m_it.finished(); }
    pixel_type dereference() const { return pixel_type(*m_ima, *m_it); }

    template <class dummy = bool>
    typename std::enable_if<has_NL::value, dummy>::type NL() const { return m_it.NL(); }


  private:
    Morpher*      m_ima;
    PixelIterator m_it;
  };


  /******************************************/
  /****             Morpher             *****/
  /******************************************/


  template <typename I>
  struct extended_by_value_image
    : morpher_base< extended_by_value_image<I>,
                    I,
                    typename std::remove_reference<I>::type::point_type,
                    typename std::remove_reference<I>::type::value_type>
  {
  private:
    typedef typename std::remove_reference<I>::type image_t;
    typedef extended_by_value_image<I>              this_t;
    friend struct morpher_core_access;

  public:
    typedef typename image_t::point_type			point_type;
    typedef typename image_t::domain_type                       domain_type;
    typedef typename image_t::value_type			value_type;
    typedef typename image_reference<image_t>::type		reference;
    typedef typename image_const_reference<image_t>::type	const_reference;

    typedef extended_by_value_image_pixel<this_t,       typename image_pixel<image_t>::type>       pixel_type;
    typedef extended_by_value_image_pixel<const this_t, typename image_const_pixel<image_t>::type> const_pixel_type;

    typedef typename image_value_range<image_t>::type		value_range;
    typedef typename image_const_value_range<image_t>::type	const_value_range;

    typedef extended_by_value_image_pixel_iterator<
      this_t, typename range_iterator<typename image_pixel_range<image_t>::type>::type > pixel_range_iterator;

    typedef extended_by_value_image_pixel_iterator<
      const this_t, typename range_iterator<typename image_const_pixel_range<image_t>::type>::type > const_pixel_range_iterator;

    typedef iterator_range<pixel_range_iterator>		pixel_range;
    typedef iterator_range<const_pixel_range_iterator>          const_pixel_range;

    typedef by_value_extension<value_type>	extension_type;


    extended_by_value_image(I&& ima, const value_type& val)
      : m_ima(std::forward<I>(ima)), m_val (val)
    {
    }

    friend
    internal::initializer<mln_concrete(I),
                          typename internal::image_init_from<extended_by_value_image>::type>
    imconcretize(const extended_by_value_image& f)
    {
      return std::move(imconcretize(f.m_ima));
    }

    template <typename V>
    friend
    internal::initializer<mln_ch_value(I, V),
                          typename internal::image_init_from<extended_by_value_image>::type>
    imchvalue(const extended_by_value_image& f)
    {
      return std::move(imchvalue<V>(f.m_ima));
    }


    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, reference>::type
    operator() (const point_type& p)
    {
      mln_precondition(m_ima.domain().has(p));
      return m_ima(p);
    }

    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, const_reference>::type
    operator() (const point_type& p) const
    {
      mln_precondition(m_ima.domain().has(p));
      return m_ima(p);
    }


    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, reference>::type
    at (const point_type& p)
    {
      if (m_ima.domain().has(p))
	return m_ima(p);
      else
	return m_val;
    }

    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, const_reference>::type
    at (const point_type& p) const
    {
      if (m_ima.domain().has(p))
        return m_ima(p);
      else
        return m_val;
    }

    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, pixel_type>::type
    pixel (const point_type& p)
    {
      // FIXME: this is wrong: it depends if we are in the domain or not
      mln_precondition(false);
      return {*this,  m_ima.pixel(p)};
    }

    template <typename = void>
    typename std::enable_if<image_traits<image_t>::accessible::value, const_pixel_type>::type
    pixel (const point_type& p) const
    {
      // FIXME: this is wrong: it depends if we are in the domain or not
      mln_precondition(false);
      return {*this,  m_ima.pixel(p)};
    }


    const_pixel_range pixels() const
    {
      return const_pixel_range(const_pixel_range_iterator(*this, m_ima.pixels().iter()));
    }

    pixel_range pixels()
    {
      return pixel_range(pixel_range_iterator(*this, m_ima.pixels().iter()));
    }

    extension_type extension() const
    {
      return extension_type (const_cast<value_type&>(m_val));
    }


  private:
    I		 m_ima;
    mln_value(I) m_val;
  };


/******************************************/
/****            Extension            ****/
/******************************************/

  template <typename V>
  struct by_value_extension
  {
    typedef V               value_type;
    typedef std::true_type  support_fill;
    typedef std::false_type support_mirror;
    typedef std::false_type support_periodize;


    by_value_extension(V& val)
    : m_val(val)
    {
    }

    void fill(const V& val)
    {
      m_val = val;
    }

  private:
    V& m_val;
  };


  template <typename I>
  extended_by_value_image<const I&>
  extend_by_value(const Image<I>& ima, const mln_value(I)& v)
  {
    return extended_by_value_image<const I&>(exact(ima), v);
  }

  template <typename I>
  extended_by_value_image<I&>
  extend_by_value(Image<I>& ima, const mln_value(I)& v)
  {
    return extended_by_value_image<I&>(exact(ima), v);
  }

  template <typename I>
  extended_by_value_image<I>
  extend_by_value(Image<I>&& ima, const mln_value(I)& v)
  {
    return extended_by_value_image<I>(move_exact(ima), v);
  }


}


#endif // ! EXTENDED_BY_VALUE_IMAGE_HPP
