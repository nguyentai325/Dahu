#ifndef MLN_CORE_IMAGE_CONSTANT_IMAGE_HPP
# define MLN_CORE_IMAGE_CONSTANT_IMAGE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/range/range.hpp>
# include <mln/core/range/iterator_range.hpp>

// FIXME: specialize for non bidirectional domain.

namespace mln
{
  template <class Domain, class V>
  struct constant_image;


  /******************************************/
  /****              Traits              ****/
  /******************************************/

  template <class Domain, class V>
  struct image_traits< constant_image<Domain, V> >
  {
    typedef std::false_type		concrete;
    typedef forward_image_tag	category;
    typedef std::true_type		accessible;
    typedef std::true_type		indexable;
    typedef mln::extension::value_extension_tag extension;
  };

  namespace internal
  {
    template <class Domain, class V>
    struct image_init_from< constant_image<Domain, V> >
    {
      typedef Domain type;
    };
  }


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <class Domain, class V>
  struct constant_image : image_base< constant_image<Domain, V>,
				      typename range_value<Domain>::type,
				      V>
  {
  private:
    typedef typename range_iterator<Domain>::type piter;

  public:
    typedef Domain				domain_type;
    typedef typename range_value<Domain>::type  point_type;
    typedef point_type				site_type;

    typedef V value_type;
    typedef V reference;
    typedef V const_reference;

    typedef unsigned size_type;
    typedef int	     difference_type;


    struct const_value_iterator : iterator_base<const_value_iterator,
						V, V>
    {
      const_value_iterator(const piter& x, const V& value)
        : m_x(x), m_value(value)
      {
      }
      void init() { m_x.init(); }
      void next() { m_x.next(); }
      bool finished() const { return m_x.finished(); }
      V dereference() const { return m_value; }
    private:
      piter m_x;
      V	    m_value;
    };

    template <class I>
    struct pixel_t : Pixel< pixel_t<I> >
    {
      template <class>
      friend struct pixel_t;

      typedef constant_image::point_type point_type;
      typedef point_type site_type;
      typedef I image_type;
      typedef V value_type;
      typedef V reference;
      typedef unsigned size_type;

      pixel_t() = default;
      pixel_t(const pixel_t&) = default;

      template <class dummy = void>
      pixel_t(const pixel_t<typename std::remove_const<I>::type>& other,
	      typename std::enable_if<std::is_const<I>::value, dummy>::type* = NULL)
	: m_ima (other.m_ima), m_p (other.m_p) , m_value (other.m_value)
      {
      }


      pixel_t(I* ima, const point_type& p, const V& val)
	: m_ima (ima), m_p (p), m_value(val)
      {
      }
      point_type  point() const { return m_p; }
      point_type  site() const { return m_p; }
      V		  val() const { return m_value; }
      I&          image() const { return *m_ima; }
      constexpr
      size_type   index() const { return 0; }
    private:
      I* m_ima;
      point_type  m_p;
      V		  m_value;
    };

    template <class I>
    struct pixel_iterator_t : iterator_base<pixel_iterator_t<I>,
					    pixel_t<I>, pixel_t<I> >
    {
    private:
      typedef I image_type;

      template <class>
      friend struct pixel_iterator_t;

    public:
      pixel_iterator_t() = default;
      pixel_iterator_t(const pixel_iterator_t&) = default;

      template <class dummy = void>
      pixel_iterator_t(const pixel_iterator_t<typename std::remove_const<I>::type>& other,
		       typename std::enable_if<std::is_const<I>::value, dummy>::type* = NULL)
	: m_ima (other.m_ima), m_x (other.m_x), m_value (other.m_value)
      {
      }

      pixel_iterator_t(I* ima, const piter& x, const V& val)
	: m_ima(ima), m_x(x), m_value(val)
      {
      }
      void	 init() { m_x.init(); }
      void	 next() { m_x.next(); }
      pixel_t<I> dereference() const { return pixel_t<I>(m_ima, *m_x, m_value); }
      bool	 finished() const { return m_x.finished(); }
    private:
      I*	      m_ima;
      piter	      m_x;
      V		      m_value;
    };



    typedef pixel_t<constant_image>		pixel_type;
    typedef pixel_t<const constant_image>	const_pixel_type;


    typedef iterator_range<const_value_iterator> value_range;
    typedef iterator_range<const_value_iterator> const_value_range;
    typedef iterator_range< pixel_iterator_t<constant_image>  >		pixel_range;
    typedef iterator_range< pixel_iterator_t<const constant_image> >	const_pixel_range;


    /************************************/
    /***  Extension                    **/
    /************************************/

    struct extension_type
    {
      typedef V               value_type;
      typedef std::false_type support_fill;
      typedef std::false_type support_mirror;
      typedef std::false_type support_periodize;
    };




    constant_image(const Domain& dom, V value)
      : m_domain(dom), m_value(value)
    {
    }

    friend
    internal::initializer<mln_concrete(constant_image),
                          Domain>
    imconcretize(const constant_image& f)
    {
      return { f.m_domain };
    }

    template <typename T>
    friend
    internal::initializer<mln_ch_value(constant_image, T),
                          Domain>
    imchvalue(const constant_image& f)
    {
      return { f.m_domain };
    }

    const domain_type&
    domain() const
    {
      return m_domain;
    }

    value_range
    values() const
    {
      return value_range( const_value_iterator(rng::iter(m_domain), m_value) );
    }

    value_range
    values()
    {
      return ((const constant_image*)this)->values();
    }

    const_pixel_range
    pixels() const
    {
      return const_pixel_range( pixel_iterator_t<const constant_image> (this, rng::iter(m_domain), m_value) );
    }

    pixel_range
    pixels()
    {
      return pixel_range( pixel_iterator_t<constant_image> (this, rng::iter(m_domain), m_value) );
    }


    V operator () (const point_type& p) const
    {
      mln_precondition(m_domain.has(p));
      (void) p;
      return m_value;
    }

    V operator () (const point_type& p)
    {
      return (*(const constant_image*)this) (p);
    }

    constexpr
    V at(const point_type&) const
    {
      return m_value;
    }

    V at(const point_type&)
    {
      return m_value;
    }

    pixel_type
    pixel(const point_type& p)
    {
      return pixel_type(this, p, m_value);
    }

    const_pixel_type
    pixel(const point_type& p) const
    {
      return const_pixel_type(this, p, m_value);
    }

    pixel_type
    pixel_at(const point_type& p)
    {
      return pixel_type(this, p, m_value);
    }

    const_pixel_type
    pixel_at(const point_type& p) const
    {
      return const_pixel_type(this, p, m_value);
    }


    constexpr
    V operator [] (size_type) const
    {
      return m_value;
    }

    V operator [] (size_type)
    {
      return m_value;
    }


    constexpr
    size_type index_of_point(const point_type&) const
    {
      return 0;
    }

    constexpr
    point_type point_at_index(size_type) const
    {
      return point_type ();
    }

    constexpr
    difference_type delta_index(const point_type&) const
    {
      return 0;
    }

    constexpr
    extension_type extension() const
    {
      return extension_type ();
    }


  private:
    Domain m_domain;
    V	   m_value;
  };

}

#endif // ! MLN_CORE_IMAGE_CONSTANT_IMAGE_HPP
