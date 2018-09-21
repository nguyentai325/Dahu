#ifndef MLN_CORE_PIXEL_UTILITY_HPP
# define MLN_CORE_PIXEL_UTILITY_HPP

# include <mln/core/image/morphers/details/morpher_core_access.hpp>
# include <mln/core/iterator/iterator_base.hpp>

namespace mln
{
  /// \brief Base class for morphed pixels thats acts like a mixin on Pix.
  template <typename Derived, typename Pix,
	    typename Morpher = typename Pix::image_type>
  struct morpher_pixel_base;


  /// \brief A pixel morpher that does nothing but changing the image binded.
  template <typename I, typename Pixel>
  struct rebinded_pixel;

  template <typename I, typename PixelIterator>
  struct rebind_pixel_iterator;

  // Meta operators for pixels
  struct get_pixel_value;
  struct get_pixel_point;

  namespace internal
  {
    /// \brief Helper class to know if a function accepts pixels or values
    template <typename I, class UnaryFunction>
    struct use_pix_helper;


    typedef mln::get_pixel_value meta_get_pixel_value;
    typedef mln::get_pixel_point meta_get_pixel_point;
  }

  /******************************************/
  /****          Implementation          ****/
  /******************************************/


  namespace internal
  {

    template <typename I, class UnaryFunction>
    struct use_pix_helper
    {
      template <typename U>
      static
      std::true_type foo(const U&, decltype( std::declval<UnaryFunction>() (std::declval<U>()) )* = NULL);

      static
      std::false_type foo(...);

      typedef decltype(foo(std::declval<mln_pixel(I)>())) type;
    };

  }

  /******************************************/
  /****          HELPER MACROS          *****/
  /******************************************/

# define MLN_PIXMORPHER_FORWARD_IF_0_(COND, F, RETURN, CV)      \
  typename std::enable_if<COND, RETURN>::type                   \
  F() CV                                                        \
  {                                                             \
    return morpher_core_access::get_pix(this).F();              \
  }

# define MLN_PIXMORPHER_FORWARD_CONST_0(F, RETURN)    MLN_PIXMORPHER_FORWARD_IF_0_(true, F, RETURN, const)


  /******************************************/
  /****        morpher pixel base        ****/
  /******************************************/

  namespace impl
  {
    // If the image of the pixel is indexable
    // we add the following typedefs/methods:
    // + typedef size_type
    // + index()
    template <typename Derived, typename Pix, typename Morpher,
	      bool is_indexable =
	      image_traits<Morpher>::indexable::value
	      >
    struct morpher_pixel_indexable;

    template <typename Derived, typename Pix, typename Morpher>
    struct morpher_pixel_indexable<Derived, Pix, Morpher, false>
    {
    };

    template <typename Derived, typename Pix, typename Morpher>
    struct morpher_pixel_indexable<Derived, Pix, Morpher, true>
    {
    private:
      friend  struct mln::morpher_core_access;
      typedef Pix         pixel_t;
      typedef Derived     derived_t;

    public:
      typedef typename Pix::size_type	size_type;

      MLN_PIXMORPHER_FORWARD_CONST_0(index, size_type);
    };

  }

  template <typename Derived, typename Pix, typename Morpher>
  struct morpher_pixel_base :
    Pixel<Derived>,
    impl::morpher_pixel_indexable<Derived, Pix, Morpher>
  {
  private:
    friend  struct morpher_core_access;
    typedef Pix         pixel_t;
    typedef Derived     derived_t;

  public:
    typedef typename Pix::value_type value_type;
    typedef typename Pix::point_type point_type;
    typedef typename Pix::site_type  site_type;
    typedef typename Pix::reference  reference;
    typedef Morpher                  image_type;

    MLN_PIXMORPHER_FORWARD_CONST_0(val, reference);
    MLN_PIXMORPHER_FORWARD_CONST_0(point, point_type);
    MLN_PIXMORPHER_FORWARD_CONST_0(site, site_type);

  protected:
    Pix&		get_morphed() { return morpher_core_access::get_pix_(this); }
    const Pix&		get_morphed() const { return morpher_core_access::get_pix_(this); }
  };


  /******************************************/
  /****          Rebinded pixel          ****/
  /******************************************/


  template <typename I, typename Pixel>
  struct rebinded_pixel :
    morpher_pixel_base< rebinded_pixel<I, Pixel>,
                        typename std::remove_reference<Pixel>::type, I>
  {
    rebinded_pixel() = default;
    rebinded_pixel(const rebinded_pixel&) = default;
    rebinded_pixel(I& ima, const Pixel& pix)
      : m_ima(&ima), m_pix(pix)
    {
    }

    template <typename J, typename Pixel2>
    rebinded_pixel(const rebinded_pixel<J, Pixel2>& other,
                   typename std::enable_if<std::is_convertible<J*, I*>::value and
                   std::is_convertible<Pixel2, Pixel>::value>::type* = NULL)
      : m_ima(other.m_ima),
	m_pix(other.m_pix)
    {
    }

    I& image() const
    {
      return *m_ima;
    }

  private:
    friend  struct morpher_core_access;

    template <typename, typename>
    friend struct rebinded_pixel;
    I*    m_ima;
    Pixel m_pix;
  };


  template <typename I, typename PixelIterator>
  struct rebind_pixel_iterator :
    iterator_base <rebind_pixel_iterator<I, PixelIterator>,
                   rebinded_pixel<I, typename PixelIterator::value_type>,
                   rebinded_pixel<I, typename PixelIterator::reference> >
  {
  private:
    typedef rebinded_pixel<I, typename PixelIterator::reference> pixel_t;

  public:
    rebind_pixel_iterator() = default;

    rebind_pixel_iterator(const rebind_pixel_iterator& other) = default;

    rebind_pixel_iterator(I& ima, const PixelIterator& pixter)
      : m_ima(&ima),
	m_pixter(pixter)
    {
    }

    template <typename J, typename PixelIterator2>
    rebind_pixel_iterator(const rebind_pixel_iterator<J, PixelIterator2>& other,
                            typename std::enable_if< std::is_convertible<J*, I*>::value and
                            std::is_convertible<PixelIterator2, PixelIterator>::value>::type* = NULL)
      : m_ima(other.m_ima),
	m_pixter(other.m_pixter)
    {
    }

    void init()
    {
      m_pixter.init();
    }

    void next()
    {
      m_pixter.next();
    }

    pixel_t dereference() const
    {
      return pixel_t(*m_ima, *m_pixter);
    }

    bool finished() const
    {
      return m_pixter.finished();
    }


  private:
    template <typename, typename>
    friend struct rebind_pixel_iterator;
    I*			m_ima;
    PixelIterator       m_pixter;
  };


  struct get_pixel_value
  {
    template <typename Pix>
    struct apply { typedef typename Pix::reference type; };

    template <typename Pix>
    typename Pix::reference
    operator() (const Pix& pix) const
    {
      return pix.val();
    }
  };

  struct get_pixel_point
  {
    template <typename Pix>
    struct apply { typedef typename Pix::point_type type; };

    template <typename Pix>
    typename Pix::point_type
    operator() (const Pix& pix) const
    {
      return pix.point();
    }
  };

} // end of namespace mln

# undef MLN_PIXMORPHER_FORWARD_0_
# undef MLN_PIXMORPHER_FORWARD_CONST_0

#endif //!MLN_CORE_PIXEL_UTILITY_HPP
