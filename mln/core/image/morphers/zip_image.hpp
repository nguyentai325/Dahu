#ifndef MLN_CORE_IMAGE_MORPHERS_ZIP_IMAGE_HPP
# define MLN_CORE_IMAGE_MORPHERS_ZIP_IMAGE_HPP

# include <tuple>
# include <mln/core/image/morphers/morpher_base.hpp>
# include <mln/core/range/zip.hpp>
# include <mln/core/range/transform.hpp>
# include <mln/core/internal/tuple_utility.hpp>


namespace mln
{

  template <class... Images>
  struct zip_image;


  template <class... Images>
  zip_image<Images...>
  imzip(Images&&... images);


  /******************************************/
  /****              Traits              ****/
  /******************************************/

  // Specilization for a single image
  template <class I_>
  struct image_traits< zip_image<I_> >
  {
  private:
    typedef typename std::remove_reference<I_>::type I;
    typedef image_traits<I> ITraits;

  public:
    typedef typename std::common_type< random_access_image_tag,
				       typename ITraits::category>::type category;

    typedef std::false_type		  concrete;
    typedef typename ITraits::accessible  accessible;
    typedef typename ITraits::indexable   indexable;
    typedef extension::none_extension_tag extension;
  };

  // Tail recursion (2 images)
  template <class I_, class J_>
  struct image_traits< zip_image<I_,J_> >
  {
  private:
    typedef typename std::remove_reference<I_>::type I;
    typedef typename std::remove_reference<J_>::type J;
    typedef image_traits<I> ITraits;
    typedef image_traits<J> JTraits;

  public:

    typedef typename std::common_type< random_access_image_tag,
				       typename ITraits::category,
				       typename JTraits::category>::type category;

    typedef std::false_type						concrete;
    typedef std::integral_constant<bool,
				   ITraits::accessible::value &&
				   JTraits::accessible::value>  accessible;
    typedef std::integral_constant<bool,
				   ITraits::indexable::value &&
				   JTraits::indexable::value>   indexable; // FIXME: maybe
    typedef extension::none_extension_tag			extension;
  };

  // Recursive definition
  template <class I_, class... Images>
  struct image_traits< zip_image<I_, Images...> >
  {
  private:
    typedef typename std::remove_reference<I_>::type I;
    typedef zip_image<Images...>		     Tail;
    typedef image_traits<I>	HTraits;
    typedef image_traits<Tail>  TTraits;

  public:
    typedef typename std::common_type< typename HTraits::category,
				       typename TTraits::category >::type category;

    typedef std::false_type					concrete;
    typedef std::integral_constant<bool,
				   HTraits::accessible::value and
				   TTraits::accessible::value>  accessible;
    typedef std::integral_constant<bool,
				   HTraits::indexable::value and
				   TTraits::indexable::value>   indexable;
    typedef extension::none_extension_tag			extension;
  };

  template <class I0, class... I>
  struct image_concrete< zip_image<I0, I...> >
  {
    typedef std::tuple<typename std::decay<mln_value(I0)>::type,
                       typename std::decay<mln_value(I)>::type...> vtype;
    typedef mln_ch_value(I0, vtype) type;
  };

  namespace internal
  {
    template <class I0, class... I>
    struct image_init_from< zip_image<I0, I...> >
    {
      typedef typename image_init_from<
        typename std::decay<I0>::type
        >::type type;
    };
  }


  /******************************************/
  /****          Implementation          ****/
  /******************************************/


  template <class... Images>
  zip_image<Images...>
  imzip(Images&&... images)
  {
    //BOOST_CONCEPT_ASSERT(Image<Images>)...;

    return zip_image<Images...>(std::forward<Images>(images)...);
  }

  namespace internal
  {

    struct meta_get_value_range
    {
      template <typename I>
      struct apply
      {
	typedef mln_vrange(I) type;
      };

      template <typename I>
      mln_vrange(I)
      operator() (I& ima) const
      {
	return ima.values();
      }
    };

    struct meta_get_const_value_range
    {
      template <typename I>
      struct apply
      {
	typedef typename image_const_value_range< typename std::remove_reference<I>::type >::type type;
      };

      template <typename I>
      typename apply<I>::type
      operator() (const I& ima) const
      {
	return ima.values();
      }
    };

    struct meta_get_pixel_range
    {
      template <typename I>
      struct apply
      {
	typedef mln_pixrange(I) type;
      };

      template <typename I>
      mln_pixrange(I)
      operator() (I& ima) const
      {
	return ima.pixels();
      }
    };

    struct meta_get_const_pixel_range
    {
      template <typename I>
      struct apply
      {
	typedef mln_cpixrange(I) type;
      };

      template <typename I>
      mln_cpixrange(I)
      operator() (const I& ima) const
      {
	return ima.pixels();
      }
    };

    template <typename P>
    struct meta_access_value
    {
      template <typename I>
      struct apply
      {
	typedef mln_reference(I) type;
      };

      template <typename I>
      mln_reference(I)
      operator() (I& ima) const
      {
	return ima(m_p);
      }

      const P& m_p;
    };

    template <typename P>
    struct meta_const_access_value
    {
      template <typename I>
      struct apply
      {
	typedef mln_creference(I) type;
      };

      template <typename I>
      mln_creference(I)
      operator() (const I& ima) const
      {
	return ima(m_p);
      }

      const P& m_p;
    };

    template <typename P>
    struct meta_access_pixel
    {
      template <typename I>
      struct apply
      {
    	typedef mln_pixel(I) type;
      };

      template <typename I>
      mln_pixel(I)
      operator() (I& ima) const
      {
    	return ima.pixel(m_p);
      }

      const P& m_p;
    };

    template <typename P>
    struct meta_const_access_pixel
    {
      template <typename I>
      struct apply
      {
    	typedef mln_cpixel(I) type;
      };

      template <typename I>
      mln_cpixel(I)
      operator() (const I& ima) const
      {
    	return ima.pixel(m_p);
      }

      const P& m_p;
    };


    template <typename P>
    struct meta_at
    {
      template <typename I>
      struct apply
      {
	typedef mln_reference(I) type;
      };

      template <typename I>
      mln_reference(I)
      operator() (I& ima) const
      {
	return ima.at(m_p);
      }

      const P& m_p;
    };

    template <typename P>
    struct meta_const_at
    {
      template <typename I>
      struct apply
      {
	typedef mln_creference(I) type;
      };

      template <typename I>
      mln_creference(I)
      operator() (const I& ima) const
      {
	return ima.at(m_p);
      }

      const P& m_p;
    };


    template <typename size_type>
    struct meta_access_index
    {
      template <typename I>
      struct apply
      {
	typedef mln_reference(I) type;
      };

      template <typename I>
      mln_reference(I)
      operator() (I& ima) const
      {
	return ima[m_i];
      }

      size_type m_i;
    };

    template <typename size_type>
    struct meta_const_access_index
    {
      template <typename I>
      struct apply
      {
	typedef mln_reference(I) type;
      };

      template <typename I>
      mln_reference(I)
      operator() (I& ima) const
      {
	return ima[m_i];
      }

      size_type m_i;
    };


    struct meta_get_pixval
    {
      template <typename Pix>
      struct apply
      {
	typedef typename std::remove_reference<Pix>::type::reference type;
      };

      template <typename Pix>
      typename Pix::reference
      operator() (const Pix& px) const
      {
	return px.val();
      }
    };

  }



  template <class... Images>
  struct zip_image : morpher_base<
    zip_image<Images...>,
    typename std::tuple_element<0, std::tuple<Images...> >::type,
    typename std::remove_reference<typename std::tuple_element<0, std::tuple<Images...> >::type>::type::point_type,
    std::tuple<typename std::decay<mln_reference(Images)>::type...>
    >
  {
    static_assert(sizeof... (Images) >= 2, "You must zip at least two images.");

  private:
    typedef zip_image<Images...> this_t;
    typedef typename std::tuple_element<0, std::tuple<Images...> >::type first_image_;
    typedef typename std::remove_reference<first_image_>::type		 I0;

    typedef std::tuple<mln_vrange(Images)...>					vrange_tuple;
    typedef std::tuple<mln_cvrange(Images)...>					cvrange_tuple;
    typedef std::tuple<mln_pixrange(Images)...>					pixrange_tuple;
    typedef std::tuple<mln_cpixrange(Images)...>				cpixrange_tuple;
    typedef std::tuple<typename range_reference<mln_pixrange(Images)>::type...>		pixel_tuple;
    typedef std::tuple<typename range_reference<mln_cpixrange(Images)>::type...>	cpixel_tuple;

    friend struct mln::morpher_core_access;

  public:
    struct pixel_type;
    struct const_pixel_type;

  private:
    struct pixel_fun_t
    {
      pixel_type
      operator() (const pixel_tuple& t) const
      {
	return pixel_type(t, m_this);
      }
      this_t* m_this;
    };

    struct const_pixel_fun_t
    {
      const_pixel_type
      operator() (const cpixel_tuple& t) const
      {
	return const_pixel_type(t, m_this);
      }
      const this_t* m_this;
    };

  public:
    typedef typename std::tuple<Images...>		images_type;
    typedef images_type					image_tuple_t;

    typedef mln_point(I0)				point_type;
    typedef typename I0::domain_type			domain_type;
    typedef std::tuple<mln_reference(Images)...>	reference;
    typedef std::tuple<mln_creference(Images)...>	const_reference;
    typedef std::tuple<typename std::decay<mln_value(Images)>::type...>            value_type;

    typedef zip_range<vrange_tuple>						value_range;
    typedef zip_range<cvrange_tuple>						const_value_range;
    typedef transformed_range<zip_range<pixrange_tuple>, pixel_fun_t>		pixel_range;
    typedef transformed_range<zip_range<cpixrange_tuple>, const_pixel_fun_t>	const_pixel_range;


    zip_image(Images&&... images)
      : m_images(std::forward_as_tuple(images...))
    {
    }

    friend
    internal::initializer<mln_concrete(zip_image),
                          typename internal::image_init_from<zip_image>::type>
    imconcretize(const zip_image& f)
    {
      using mln::imchvalue;
      return std::move(imchvalue<value_type>(std::get<0>(f.m_images)));
    }

    template <typename V>
    friend
    internal::initializer<mln_ch_value(zip_image, V),
                          typename internal::image_init_from<zip_image>::type>
    imchvalue(const zip_image& f)
    {
      using mln::imchvalue;
      return std::move(imchvalue<V>(std::get<0>(f.m_images)));
    }

    images_type& images()
    {
      return m_images;
    }

    const images_type& images() const
    {
      return m_images;
    }


    const_value_range
    values() const
    {
      auto t = internal::tuple_transform(m_images, internal::meta_get_const_value_range ());
      return const_value_range(t);
    }

    value_range
    values()
    {
      auto t = internal::tuple_transform(m_images, internal::meta_get_value_range ());
      return value_range(t);
    }

    const_pixel_range
    pixels() const
    {
      auto t = internal::tuple_transform(m_images, internal::meta_get_const_pixel_range ());
      const_pixel_fun_t fun {this};
      return const_pixel_range(t, fun);
    }

    pixel_range
    pixels()
    {
      auto t = internal::tuple_transform(m_images, internal::meta_get_pixel_range ());
      pixel_fun_t fun {this};
      return pixel_range(t, fun);
    }

    reference
    operator() (const point_type& p)
    {
      return internal::tuple_transform(m_images, internal::meta_access_value<point_type> {p} );
    }

    const_reference
    operator() (const point_type& p) const
    {
      return internal::tuple_transform(m_images, internal::meta_const_access_value<point_type> {p});
    }

    reference
    at (const point_type& p)
    {
      return internal::tuple_transform(m_images, internal::meta_at<point_type> {p} );
    }

    const_reference
    at (const point_type& p) const
    {
      return internal::tuple_transform(m_images, internal::meta_const_at<point_type> {p});
    }

    pixel_type
    pixel (const point_type& p)
    {
      auto t = internal::tuple_transform(m_images, internal::meta_access_pixel<point_type> {p} );
      return pixel_type(t, this);
    }

    const_pixel_type
    pixel (const point_type& p) const
    {
      auto t = internal::tuple_transform(m_images, internal::meta_const_access_pixel<point_type> {p} );
      return const_pixel_type(t, this);
    }

    template <typename image_t = I0>
    reference
    operator[] (typename image_t::size_type i)
    {
      return internal::tuple_transform(m_images, internal::meta_access_index<typename image_t::size_type> {i});
    }

    template <typename image_t = I0>
    const_reference
    operator[] (typename image_t::size_type i) const
    {
      return internal::tuple_transform(m_images, internal::meta_const_access_index<typename image_t::size_type> {i});
    }

  private:
    const I0&
    get_morphed() const
    {
      return std::get<0>(m_images);
    }

    I0&
    get_morphed()
    {
      return std::get<0>(m_images);
    }

    std::tuple<Images...>	m_images;
  };

  template <class... Images>
  struct zip_image<Images...>::pixel_type
    : morpher_pixel_base< pixel_type,
			  typename std::remove_reference<
			    typename std::tuple_element<0, pixel_tuple>::type>::type >
  {
  public:
    typedef zip_image<Images...>			image_type;

  private:
    typedef pixel_tuple					pixel_tuple_t;
    typedef typename std::remove_reference<
      typename std::tuple_element<0, pixel_tuple_t>::type>::type first_pixel_t;
    friend struct const_pixel_type;

  public:
    typedef std::tuple<mln_reference(Images)...>	reference;
    typedef reference					value_type;
    typedef typename image_type::point_type		point_type;
    typedef typename image_type::site_type		site_type;

    pixel_type(const pixel_tuple_t& px, image_type* ima)
      : m_pix_tuple(px), m_ima(ima)
    {
    }

    reference val() const
    {
      auto x = internal::tuple_transform(m_pix_tuple, internal::meta_get_pixval());
      return x;
    }

    point_type point() const
    {
      return std::get<0>(m_pix_tuple).point();
    }

    point_type site() const
    {
      return std::get<0>(m_pix_tuple).site();
    }

    image_type& image() const
    {
      return *m_ima;
    }

    template <typename image_t = image_type>
    typename std::enable_if<
      image_traits<image_t>::indexable::value,
      typename image_t::size_type>::type
    index() const
      {
	// FIXME: add assertion to check dynamically if it really indexable
	return std::get<0>(m_pix_tuple).index();
      }

  private:
    const first_pixel_t& get_morphed() const { return std::get<0>(m_pix_tuple); };
    first_pixel_t& get_morphed() { return std::get<0>(m_pix_tuple); };

    pixel_tuple_t m_pix_tuple;
    image_type*   m_ima;
  };


  template <class... Images>
  struct zip_image<Images...>::const_pixel_type
    : morpher_pixel_base< const_pixel_type,
			  typename std::remove_reference<
			    typename std::tuple_element<0, cpixel_tuple>::type>::type >
  {
  public:
    typedef const zip_image<Images...>			image_type;

  private:
    typedef cpixel_tuple				pixel_tuple_t;
    typedef typename std::remove_reference<
      typename std::tuple_element<0, pixel_tuple_t>::type>::type first_pixel_t;

  public:
    typedef std::tuple<mln_creference(Images)...>	reference;
    typedef reference					value_type;
    typedef typename image_type::point_type		point_type;
    typedef typename image_type::site_type		site_type;


    const_pixel_type(const pixel_tuple_t& px, image_type* ima)
      : m_pix_tuple(px), m_ima(ima)
    {
    }


    // interop
    const_pixel_type(const pixel_type& other)
      : m_pix_tuple(other.m_pix_tuple), m_ima(other.m_ima)
    {
    }

    reference val() const
    {
      auto x = internal::tuple_transform(m_pix_tuple, internal::meta_get_pixval());
      return x;
    }

    point_type point() const
    {
      return std::get<0>(m_pix_tuple).point();
    }

    point_type site() const
    {
      return std::get<0>(m_pix_tuple).site();
    }

    image_type& image() const
    {
      return *m_ima;
    }

    template <typename image_t = image_type>
    typename std::enable_if<
      image_traits<image_t>::indexable::value,
      typename image_t::size_type>::type
    index() const
      {
	// FIXME: add assertion to check dynamically if it really indexable
	return std::get<0>(m_pix_tuple).index();
      }

  private:
  private:
    const first_pixel_t& get_morphed() const { return std::get<0>(m_pix_tuple); };
    first_pixel_t& get_morphed() { return std::get<0>(m_pix_tuple); };

    pixel_tuple_t	m_pix_tuple;
    image_type*		m_ima;
  };


}

#endif // ! MLN_CORE_IMAGE_MORPHERS_ZIP_IMAGE_HPP
