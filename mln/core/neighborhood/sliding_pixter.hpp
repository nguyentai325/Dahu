#ifndef MLN_CORE_NEIGHBORHOOD_SLIDING_PIXTER_HPP
# define MLN_CORE_NEIGHBORHOOD_SLIDING_PIXTER_HPP

# include <mln/core/range/range.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/iterator/iterator_base.hpp>

namespace mln
{


  /// \brief Define a pixel iterator over a siteset centered on a pixel.
  ///
  /// Define an iterator over a siteset centered on a pixel.
  ///
  /// \p Pixel can be a pointer type, in that case, the iterator will
  /// be binded to the Pixel.
  /// \p Pixel can be an itetaror, in that case, the current iterator will
  /// be binded to the iterator.
  ///
  /// Note that the siteset is not copied and
  /// its lifetime should be longer that the iterator's one.
  ///
  /// Internally, the iterator is optimized so that:
  /// * if I is raw, it uses pointers
  /// * if I is indexable, it uses indexes
  /// * otherwise, it uses points.
  ///
  /// \tparam Point
  /// \tparam SiteSet
  /// \pre The pixel's image must be accessible or indexable.
  template <class Pixel, class SiteSet>
  struct sliding_pixter;




  /******************************************/
  /****          Facade                  ****/
  /******************************************/


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {

    template <typename Px>
    using image_t = typename mln::iterator_traits<Px>::value_type::image_type;

    // Forward:
    // \tparam Pixel must be a pointer to a pixel or a pixel iterator.
    template <class Pixel,
	      class SiteSet,
	      class Image = image_t<Pixel>,
	      class Enable = void>
    struct sliding_pixter_base;



  /******************************************/
  /****          Pixel Types          ****/
  /******************************************/


  /// Sliding pixel through point access
  /// Require the image to be accessible
  /// \tparam Px is either a pointer to pixel or pixel iterator
  template <class Px, class SiteSet>
  struct sliding_pixter_pixel_site :
    Pixel< sliding_pixter_pixel_site<Px, SiteSet> >
  {
  private:
    typedef typename mln::iterator_traits<Px>::value_type T;

  public:
    friend struct sliding_pixter_base<Px, SiteSet>;

    typedef typename T::point_type	point_type;
    typedef typename T::site_type	site_type;
    typedef typename T::value_type	value_type;
    typedef typename T::reference	reference;
    typedef typename T::image_type	image_type;

    sliding_pixter_pixel_site() = default;
    sliding_pixter_pixel_site(const Px& pix, const SiteSet& s)
      : m_pix(pix), m_it(rng::iter(s))
    {
    }

    point_type  point() const { return m_pix->point() + *m_it; }
    point_type  site()  const { return m_pix->point() + *m_it; }
    reference   val() const { return m_pix->image().at (this->site()); }
    image_type& image() const { return m_pix->image(); }


  private:
    typename std::conditional<std::is_pointer<Px>::value,
			      const Px, const Px&>::type m_pix;
    typename range_const_iterator<SiteSet>::type m_it;
  };


  /// Sliding pixel through indexes
  /// Require the image to be indexable
  template <class Px, class SiteSet>
  struct sliding_pixter_pixel_index :
    Pixel< sliding_pixter_pixel_index<Px, SiteSet> >
  {
  private:
    typedef typename mln::iterator_traits<Px>::value_type T;
    typedef typename T::image_type::difference_type difference_type;

  public:
    friend struct sliding_pixter_base<Px, SiteSet>;

    typedef typename T::point_type	point_type;
    typedef typename T::site_type	site_type;
    typedef typename T::value_type	value_type;
    typedef typename T::reference	reference;
    typedef typename T::image_type	image_type;
    typedef typename T::size_type	size_type;

    sliding_pixter_pixel_index() = default;
    sliding_pixter_pixel_index(const Px& pix,
			       const SiteSet& pset,
			       const difference_type* indexes)
      : m_pix(pix),
	m_site_it(rng::iter(pset)),
	m_index_it(indexes)
    {
    }

    point_type  point() const { return m_pix->point() + *m_site_it; }
    point_type  site()  const { return m_pix->point() + *m_site_it; }
    image_type& image() const { return m_pix->image(); }
    size_type	  index() const { return m_pix->index() + *m_index_it; }
    reference   val() const   { return m_pix->image()[this->index()]; }

  private:
    typename std::conditional<std::is_pointer<Px>::value,
			      const Px, const Px&>::type m_pix;
    typename range_const_iterator<SiteSet>::type  m_site_it;
    const difference_type* m_index_it;
  };


  /// Sliding pixel through offsets
  /// Require the image to be a raw_image
  template <class Px, class SiteSet>
  struct sliding_pixter_pixel_pointer :
    Pixel< sliding_pixter_pixel_pointer<Px, SiteSet> >
  {
  private:
    typedef typename mln::iterator_traits<Px>::value_type T;
    typedef typename T::image_type::difference_type difference_type;

  public:
    friend struct sliding_pixter_base<Px, SiteSet>;

    typedef typename T::point_type	point_type;
    typedef typename T::site_type	site_type;
    typedef typename T::value_type	value_type;
    typedef typename T::reference	reference;
    typedef typename T::image_type	image_type;
    typedef typename T::size_type	size_type;

    sliding_pixter_pixel_pointer() = default;
    sliding_pixter_pixel_pointer(const Px& pix,
				 const SiteSet& pset)
      : m_pix(pix),
	m_site_it(rng::iter(pset))
    {
    }

    point_type    point() const { return m_pix->point() + *m_site_it; }
    point_type    site()  const { return m_pix->point() + *m_site_it; }
    image_type&   image() const { return m_pix->image(); }
    size_type	  index() const { return m_pix->index() + m_delta_index; }
    reference     val() const
    {
      typedef typename std::remove_reference<reference>::type T;
      typedef typename std::conditional<
	std::is_const<T>::value, const char*, char*>::type buffer_ptr_t;
      typedef typename std::add_pointer<T>::type ptr_t;

      return * reinterpret_cast<ptr_t>(reinterpret_cast<buffer_ptr_t>(&(m_pix->val()))
				       + m_delta_offset);
    }

  private:
    typename std::conditional<std::is_pointer<Px>::value,
			      const Px, const Px&>::type m_pix;
    typename range_const_iterator<SiteSet>::type   m_site_it;
    difference_type m_delta_index;
    difference_type m_delta_offset;
    //const difference_type* m_index_it;
    //const difference_type* m_offset_it;
  };



    /******************************************/
    /****          Pixter Types          ****/
    /******************************************/


    //
    // Specialization for accessible only images.
    template <class Pixel, class SiteSet, class Image>
    struct sliding_pixter_base<
      Pixel,
      SiteSet,
      Image,
      typename std::enable_if< image_traits<Image>::accessible::value and
			       !image_traits<Image>::indexable::value
			       >::type>
    : iterator_base< sliding_pixter<Pixel, SiteSet>,
		     sliding_pixter_pixel_site<Pixel, SiteSet>,
		     const sliding_pixter_pixel_site<Pixel, SiteSet>&
		     >
    {

      sliding_pixter_base() = default;
      sliding_pixter_base(const Pixel& p, const SiteSet& s)
	: m_pix (p, s)
      {
      }

      void init() { m_pix.m_it.init(); }
      void next() { m_pix.m_it.next(); }
      bool finished() const { return m_pix.m_it.finished(); }

      const sliding_pixter_pixel_site<Pixel, SiteSet>&
      dereference() const { return m_pix; }

    private:
      sliding_pixter_pixel_site<Pixel, SiteSet> m_pix;
    };

    /***************************************************************/
    /****          Pixter Types: Spe for indexable images       ****/
    /***************************************************************/

    // Specialization for indexable images not raw.
    // with a static siteset
    template <class Pixel, class Image, std::size_t N, typename P>
    struct sliding_pixter_base
    < Pixel,
      std::array<P, N>,
      Image,
      typename std::enable_if< image_traits<Image>::indexable::value and
			       !std::is_same<typename image_traits<Image>::category,
					     raw_image_tag>::value>::type >
    : iterator_base< sliding_pixter<Pixel, std::array<P, N> >,
		     sliding_pixter_pixel_index< Pixel, std::array<P, N> >,
		     const sliding_pixter_pixel_index< Pixel, std::array<P, N> >&
		     >
    {
    private:
      typedef std::array<P, N> S;
      typedef std::array<typename Image::difference_type, N> C;

    public:
      sliding_pixter_base() = default;
      sliding_pixter_base(const Pixel& px, const S& s)
	: m_pix (px, s, &m_index_set[0])
      {
	Image& ima = px->image();
	auto it = rng::iter(s);
	it.init();
	for (unsigned i = 0; i < N; ++i, it.next())
	  m_index_set[i] = ima.delta_index(*it);
      }

      void init()
      {
	m_i = 0;
	m_pix.m_index_it = &m_index_set[0];
	m_pix.m_site_it.init();
      }

      void next()
      {
	++m_i;
	++m_pix.m_index_it;
	m_pix.m_site_it.next();
      }

      bool finished() const
      {
	return m_i >= N;
      }

      const sliding_pixter_pixel_index<Pixel, S>&
      dereference() const
      {
	return m_pix;
      }

    private:
      unsigned m_i;
      C m_index_set;
      sliding_pixter_pixel_index<Pixel, S> m_pix;
    };

    // Specialization for indexable images not raw.
    // with a dynamic siteset
    template <class Pixel, class SiteSet, class Image>
    struct sliding_pixter_base
    < Pixel,
      SiteSet,
      Image,
      typename std::enable_if< image_traits<Image>::indexable::value and
			       !std::is_same<typename image_traits<Image>::category,
					     raw_image_tag>::value>::type >
    : iterator_base< sliding_pixter<Pixel, SiteSet>,
		     sliding_pixter_pixel_index<Pixel, SiteSet>,
		     const sliding_pixter_pixel_index<Pixel, SiteSet>&
		     >
    {
      sliding_pixter_base() = default;

      sliding_pixter_base(const Pixel& px, const SiteSet& s)
	: m_index_set(rng::size(s)),
	  m_pix (px, s, &m_index_set[0])
      {
	Image& ima = px->image();
	auto it = rng::iter(s);
	it.init();
	unsigned n = m_index_set.size();
	for (unsigned i = 0; i < n; ++i, it.next())
	  m_index_set[i] = ima.delta_index(*it);
      }

      void init()
      {
	m_pix.m_index_it = &m_index_set[0];
	m_pix.m_site_it.init();
      }

      void next()
      {
	++m_pix.m_index_it;
	m_pix.m_site_it.next();
      }

      bool finished() const
      {
	return m_pix.m_index_it == &m_index_set[0] + m_index_set.size();
      }

      const sliding_pixter_pixel_index<Pixel, SiteSet>&
      dereference() const
      {
	return m_pix;
      }

    private:
      typedef std::vector<typename Image::difference_type> C;
      C m_index_set;
      sliding_pixter_pixel_index<Pixel, SiteSet> m_pix;
    };

    /***************************************************************/
    /****          Pixter Types: Spe for Raw images             ****/
    /***************************************************************/

    // Specialization for raw images.
    // with a static siteset
    template <class Pixel, class Image, std::size_t N, typename P>
    struct sliding_pixter_base
    < Pixel, std::array<P, N>, Image,
      typename std::enable_if< std::is_same<typename image_traits<Image>::category,
					    raw_image_tag>::value>::type >
    : iterator_base< sliding_pixter<Pixel, std::array<P, N> >,
		     sliding_pixter_pixel_pointer<Pixel, std::array<P, N> >,
		     const sliding_pixter_pixel_pointer<Pixel, std::array<P, N> >&
		     >
    {
    private:
      typedef std::array<P, N> S;
      typedef std::array<typename Image::difference_type, N> I;
      typedef std::array<typename Image::difference_type, N> O;

    public:
      sliding_pixter_base() = default;

      sliding_pixter_base(const Pixel& px, const S& s)
	: m_pix (px, s)
      {
	Image& ima = px->image();
	auto it = rng::iter(s);
	it.init();
	for (unsigned i = 0; i < N; ++i, it.next()) {
	  m_index_set[i] = ima.delta_index(*it);
	  m_offset_set[i] = ima.delta_offset(*it);
	}

      }

      void init()
      {
	m_i = 0;
	m_pix.m_site_it.init();
	if (0 < N)
	  {
	    m_pix.m_delta_index  = m_index_set[0];
	    m_pix.m_delta_offset = m_offset_set[0];
	  }
      }

      void next()
      {
	++m_i;
	m_pix.m_site_it.next();
	if (m_i < N)
	  {
	    m_pix.m_delta_index = m_index_set[m_i];
	    m_pix.m_delta_offset = m_offset_set[m_i];
	  }
      }

      bool finished() const
      {
	return m_i >= N;
      }

      const sliding_pixter_pixel_pointer<Pixel, S>&
      dereference() const
      {
	return m_pix;
      }

    private:
      unsigned m_i;
      I m_index_set;
      O m_offset_set;
      sliding_pixter_pixel_pointer<Pixel, S> m_pix;
    };


    // Specialization for raw images.
    // with a static siteset
    template <class Pixel, class SiteSet, class Image>
    struct sliding_pixter_base
    < Pixel, SiteSet, Image,
      typename std::enable_if< std::is_same<typename image_traits<Image>::category,
					    raw_image_tag>::value>::type >
    : iterator_base< sliding_pixter<Pixel, SiteSet>,
		     sliding_pixter_pixel_pointer<Pixel, SiteSet>,
		     const sliding_pixter_pixel_pointer<Pixel, SiteSet>&
		     >
    {
    private:
      typedef std::vector<typename Image::difference_type> I;
      typedef std::vector<typename Image::difference_type> O;

    public:
      sliding_pixter_base() = default;

      sliding_pixter_base(const Pixel& px, const SiteSet& s)
	: m_sz(rng::size(s)),
	  m_index_set(m_sz),
	  m_offset_set(m_sz),
	  m_pix (px, s)
      {
	Image& ima = px->image();
	auto it = rng::iter(s);
	it.init();
	for (unsigned i = 0; i < m_sz; ++i, it.next()) {
	  m_index_set[i] = ima.delta_index(*it);
	  m_offset_set[i] = ima.delta_offset(*it);
	}
      }

      void init()
      {
	m_i = 0;
	m_pix.m_site_it.init();
	if (0 < m_sz)
	  {
	    m_pix.m_delta_index = m_index_set[0];
	    m_pix.m_delta_offset = m_offset_set[0];
	  }
      }

      void next()
      {
	++m_i;
	m_pix.m_site_it.next();
	if (m_i < m_sz)
	  {
	    m_pix.m_delta_index  = m_index_set[m_i];
	    m_pix.m_delta_offset = m_offset_set[m_i];
	  }
      }

      bool finished() const
      {
	return m_i >= m_sz;
      }

      const sliding_pixter_pixel_pointer<Pixel, SiteSet>&
      dereference() const
      {
	return m_pix;
      }

    private:
      const unsigned m_sz;
      unsigned m_i;
      I m_index_set;
      O m_offset_set;
      sliding_pixter_pixel_pointer<Pixel, SiteSet> m_pix;
    };

  }

  template <class Pixel, class SiteSet>
  struct sliding_pixter : internal::sliding_pixter_base<Pixel, SiteSet>
  {
  private:
    typedef typename mln::iterator_traits<Pixel>::value_type::image_type I;
    static_assert(image_traits<I>::accessible::value or
		  image_traits<I>::indexable::value,
		  "You cannot set a neighborhood on a pixel whose image is not"
		  "accessible <ima(p)> nor indexable <ima[i]>");

  public:
    sliding_pixter(const Pixel& px, const SiteSet& s)
      : internal::sliding_pixter_base<Pixel, SiteSet>(px, s)
    {
    }
  };

}

#endif // ! MLN_CORE_NEIGHBORHOOD_SLIDING_PIXTER_HPP
