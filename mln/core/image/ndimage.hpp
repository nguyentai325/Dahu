#ifndef MLN_CORE_IMAGE_NDIMAGE_HH
# define MLN_CORE_IMAGE_NDIMAGE_HH

/// \file

# include <mln/core/image/image.hpp>
# include <mln/core/domain/box.hpp>
# include <mln/core/memory.hpp>
# include <mln/core/assert.hpp>
# include <mln/core/image_traits.hpp>
# include <mln/core/image_category.hpp>
# include <mln/core/image/ndimage_iter.hpp>



namespace mln
{

  /// \brief base class for ndimage
  template <typename T, unsigned dim, typename E> struct ndimage_base;

  // FWD
  //template <typename I, typename T> struct ndimage_iter;
  //template <typename I, typename T> struct ndimage_rev_iter;
  //template <typename T, unsigned dim, typename I> struct ndimage_pixel_iterator;
  //template <typename T, unsigned dim, typename I> struct ndimage_rev_pixel_iterator;
  template <typename T, unsigned dim, typename I> struct ndimage_pixel;
  template <typename T, unsigned dim, typename E> struct ndimage_base;


  /******************************************/
  /****              Traits              ****/
  /******************************************/


  template <typename T, unsigned dim, typename E>
  struct image_traits< ndimage_base<T, dim, E> >
  {
    typedef raw_image_tag               category;
    typedef std::true_type              accessible;
    typedef std::true_type		indexable;
    typedef std::true_type		concrete;
    typedef std::true_type		shallow_copy;
    typedef std::true_type		has_border;
    typedef mln::extension::border_extension_tag        extension;
  };

  /******************************************/
  /****            Definition            ****/
  /******************************************/

  namespace internal
  {
    template <typename T, unsigned dim>
    struct ndimage_data
    {
      // Constructor for library managed memory
      ndimage_data(size_t* shape_, unsigned border, T v = T());
      ~ndimage_data();

      // These elements are used to ensure a correct destruction
      size_t shape[dim];
      size_t strides[dim];

      // These elements are used to free memory
      size_t nbytes;
      char*  buffer;

    private:
      ndimage_data(const ndimage_data&);
    };

    template <typename T, unsigned dim>
    struct ndimage_extension
    {
    private:
      typedef point<short, dim> P;

    public:
      typedef T value_type;
      typedef std::true_type    support_fill;
      typedef std::true_type   support_mirror;
      typedef std::false_type   support_periodize;


      ndimage_extension(char* ptr, const std::size_t* strides, const P& shp, int border);
      void fill(const T& v);
      void mirror();

    private:
      template <unsigned d>
      typename std::enable_if<(d < dim)>::type
      _fill(char* ptr, const T& v);

      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      _fill(char* ptr, const T& v);

      template <unsigned d>
      typename std::enable_if<(d < dim)>::type
      _fillall(char* ptr, const T& v);

      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      _fillall(char* ptr, const T& v);

      template <unsigned d>
      typename std::enable_if<(d < dim)>::type
      _mirror(char* ptr);

      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      _mirror(char* ptr);


      template <unsigned d>
      typename std::enable_if<(d < dim)>::type
      _copy_line(char* src, char* dest);

      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      _copy_line(char* src, char* dest);


    private:
      size_t              m_strides[dim];
      P                   m_shp;
      char*               m_ptr;
      int                 m_border;
    };


  } // end of namespace mln::internal


  template <typename T, unsigned dim, typename E>
  struct ndimage_base
#ifndef MLN_DOXYGEN
    : image_base<E, point<short, dim>, T,
                 ndimage_pixel<T, dim, E>,
                 ndimage_pixel<const T, dim, const E> >
#endif
  {
  private:
    typedef ndimage_base<T, dim, E>                             this_type;

    template <class, unsigned, class>
    friend struct ndimage_base;

    template <unsigned d, typename T1, typename E1, typename T2, typename E2>
    friend
    bool are_indexes_compatible(const ndimage_base<T1, d, E1>& self,
                                const ndimage_base<T2, d, E2>& other);

  public:
    /// \name Image point/value/pixel types
    /// \{

    /// \copydoc image::site_type
    typedef point<short,dim>                                    site_type;

    /// \copydoc image::point_type
    typedef point<short,dim>                                    point_type;

    /// \copydoc image::pixel_type
    typedef ndimage_pixel<T, dim, E>                            pixel_type;

    /// \copydoc image::const_pixel_type
    typedef ndimage_pixel<const T, dim, const E>                const_pixel_type;

    /// \copydoc image::value_type
    typedef T                   value_type;

    /// \copydoc image::reference
    typedef T&                  reference;

    /// \copydoc image::const_reference
    typedef const T&            const_reference;

    /// \copydoc image::difference_type
    typedef int                 difference_type;

    /// \copydoc image::size_type
    typedef unsigned            size_type;
    typedef unsigned            index_type;

    typedef T*                  pointer;
    typedef const T*            const_pointer;
    /// \}


    /// \name Image Ranges Types
    /// \{

    /// \copydoc image::domain_type
    typedef box<short, dim>                                             domain_type;

    /// \copydoc image::value_range
    typedef ndimage_value_range<this_type, T>                           value_range;

    /// \copydoc image::const_value_range
    typedef ndimage_value_range<const this_type, const T>               const_value_range;

    /// \copydoc image::pixel_range
    typedef ndimage_pixel_range<this_type, T>                           pixel_range;

    /// \copydoc image::const_pixel_range
    typedef ndimage_pixel_range<const this_type, const T>               const_pixel_range;
    /// \}

    enum { ndim = dim};

    // As an Image
    // Constructors
    // \{
    explicit ndimage_base(unsigned border = 3);
    explicit ndimage_base(const domain_type& domain, unsigned border = 3, T v = T());

    template <typename U, class E2>
    ndimage_base(const ndimage_base<U, dim, E2>& other, mln::init);

    template <typename U, class E2>
    ndimage_base(const ndimage_base<U, dim, E2>& other, unsigned border, T v = T());
    // \}

    /// \brief Constructors from external sources
    /// \{
    static
    E
    from_buffer(void* buffer, const domain_type& domain, bool copy = false);

    static
    E
    from_buffer(void* buffer, const domain_type& domain, const size_t* strides, bool copy = false);
    /// \}

    /// \name Accession Operators
    /// \{

    /// \copydoc image::operator()(const site_type& p) const
    reference		operator() (const site_type& p);

    /// \copydoc image::operator()(const site_type& p) const
    const_reference	operator() (const site_type& p) const;

    /// \copydoc image::operator[](size_type i) const
    reference		operator[] (size_type i);

    /// \copydoc image::operator[](size_type i) const
    const_reference	operator[] (size_type i) const;

    /// \copydoc image::at(const site_type& p) const
    reference		at (const site_type& p);

    /// \copydoc image::at(const site_type& p) const
    const_reference	at (const site_type& p) const;

    /// \}

    /// \name Pixel utilities
    /// \{

    /// \copydoc image::pixel_at(const point_type&) const
    pixel_type          pixel_at(const point_type& p);

    /// \copydoc image::pixel_at(const point_type&) const
    const_pixel_type    pixel_at(const point_type& p) const;

    /// \copydoc image::pixel(const point_type&) const
    pixel_type          pixel(const point_type& p);

    /// \copydoc image::pixel(const point_type&) const
    const_pixel_type    pixel(const point_type& p) const;

    /// \}


    typedef typename value_range::iterator                      value_iterator;
    typedef typename value_range::reverse_iterator              reverse_value_iterator;
    typedef typename const_value_range::iterator                const_value_iterator;
    typedef typename const_value_range::reverse_iterator		const_reverse_value_iterator;
    typedef typename pixel_range::iterator                      pixel_iterator;
    typedef typename pixel_range::reverse_iterator              reverse_pixel_iterator;
    typedef typename const_pixel_range::iterator                const_pixel_iterator;
    typedef typename const_pixel_range::reverse_iterator		const_reverse_pixel_iterator;

    /// \name Image Ranges
    /// \{

    /// \copydoc image::domain()
    const domain_type&          domain() const;

    /// \copydoc image::values()
    value_range                 values();

    /// \copydoc image::values() const
    const_value_range           values() const;

    /// \copydoc image::pixels()
    pixel_range                 pixels();

    /// \copydoc image::pixels() const
    const_pixel_range           pixels() const;
    /// \}

    /// \name Concrete-related Image Methods
    /// \{

    /// \brief Resize the image to fit \p domain.
    ///
    /// Resize the image w.r.t to \p domain and \p border. The image is
    /// initialized with value  \p v.
    /// \warning The old values are not kept.
    /// \param domain The new domain of the image.
    /// \param border The new border of the image.
    /// \param v The initialization value of the image. If \p is given, values
    /// are copy constructed from \p v, otherwise, they are default-initialized.
    void resize(const domain_type& domain, unsigned border = 3, T v = T());

    /// \brief Re-indexation of the image.
    ///
    /// Reindex the image so that the index of the first point is \p
    /// index_first.
    /// \param index_first
    /// \postcondition `this->index_of_point(domain().pmin) == index_first`
    void reindex(size_type index_first)
    {
      std::ptrdiff_t diff = index_first - m_index_first;
      m_ptr_origin -= diff;
      m_index_first += diff;
      m_index_last += diff;
    }
    /// \}




    /// \name Index-related methods
    /// \{

    /// \copydoc image::index_of_point(const point_type&) const
    size_type       index_of_point(const point_type& p) const;

    /// \copydoc image::point_at_index(size_type i) const
    point_type		point_at_index(size_type i) const;

    /// \copydoc image::delta_index(const point_type&) const
    difference_type delta_index(const point_type& p) const;

    /// \}



    // As a Raw Image
    const std::size_t*       strides() const;
    int border() const { return border_; }


    // Specialized algorithm
    template <typename T_, unsigned dim_, typename E_, typename Domain_>
    friend typename std::enable_if<std::is_convertible<Domain_, typename ndimage_base<T_, dim_, E_>::domain_type>::value, E_>::type
    make_subimage(ndimage_base<T_, dim_, E_>&, const Domain_& domain);
    // template <typename T_, unsigned dim_, typename E_, typename Domain_>
    // friend E_ make_subimage(ndimage_base<T_, dim_, E_>&, const Domain_& domain);
    // template <typename T_, unsigned dim_, typename E_, typename Domain_>
    // friend typename E_ make_subimage(ndimage_base<T_, dim_, E_>&&, const Domain_& domain);


    // As an Indexable Image
    const size_t*	index_strides() const     { return &m_index_strides[0]; }

    // template <typename O>
    // bool friend_index_compatible(ndimage_base self, const Image<O>& other) const;






    difference_type delta_offset(const point_type& p) const
    {
      difference_type idx = 0;
      for (unsigned i = 0; i < dim; ++i)
        idx += p[i] * strides_[i];
      return idx;
    }


    // Extension
    typedef internal::ndimage_extension<T, dim> extension_type;
    extension_type extension() const;

  private:
    static
    E
    from_buffer_copy_(void* buffer, const domain_type& domain, const size_t* strides);

    static
    E
    from_buffer_extern_(void* buffer, const domain_type& domain, const size_t* strides);




  protected:
    friend struct ndimage_pixel<T, dim, E>;
    friend struct ndimage_pixel<const T, dim, const E>;
    template <typename, typename, unsigned, typename> friend struct indexible_ndimage_base;
    template <typename, typename> friend struct ndimage_value_range;
    template <typename, typename> friend struct ndimage_pixel_range;

    void resize_(const domain_type& domain, unsigned border = 3, T v = T());

    domain_type	domain_;	///< Domain of image
    //#ifndef MLN_NDEBUG
    domain_type vbox_;
    //#endif

    std::array<size_t, dim>	strides_;	///< Strides in bytes
    std::shared_ptr< internal::ndimage_data<T, dim> > data_;
    int		border_;
    char*	ptr_;           ///< Pointer to the first element in the domain
    char*	last_;          ///< Pointer to the last element in the domain (not past-the-end)

    T*					m_ptr_origin;		///< Pointer to the first element
    std::array<std::size_t, dim>	m_index_strides;	///< Strides in number of elements (including the border)
    size_t				m_index_first;          ///< index of pmin
    size_t				m_index_last;           ///< index of pmax
  };

  /******************************/
  /** Free function Impl        */
  /******************************/


  template <typename T>
  struct image3d : ndimage_base<T, 3, image3d<T> >
  {
  protected:
    typedef ndimage_base<T, 3, image3d<T> > base;
    typedef typename base::domain_type domain_type;

  public:
    explicit image3d (unsigned border = 3) : base (border) {}

    explicit image3d(const domain_type& domain, unsigned border = 3)
      : base(domain, border)
    {
    }

    image3d(short nslices, short nrows, short ncols, unsigned border = 3)
      : base( (box<short,3>) {{0,0,0},{nslices, nrows, ncols}}, border)
    {
    }
  };


  namespace internal
  {
    template <typename T, unsigned dim>
    ndimage_data<T, dim>::ndimage_data(size_t* shape_, unsigned border, T v)
    {
      for (unsigned i = 0; i < dim; ++i)
        shape[i] = shape_[i] + 2 * border;

      strides[dim-1] = sizeof(T);

      // Each row / page ... are 16 bytes aligned
      unsigned ndim = dim;

#ifdef MLN_128B_ALIGNMENT
      if (ndim >= 2)
        {
          strides[dim-2] = ((shape[dim-1] * sizeof(T)) & ~(size_t)15) + (size_t) 16;
          for (int i = dim-3; i >= 0; --i)
            strides[i] = shape[i+1] * strides[i+1];
        }
#else
      if (ndim >= 2)
        {
          strides[dim-2] = shape[dim-1] * sizeof(T);
          for (int i = dim-3; i >= 0; --i)
            strides[i] = shape[i+1] * strides[i+1];
        }
#endif

      // Allocate data
      nbytes = shape[0] * strides[0];
      buffer = (char*) mln::aligned_malloc(nbytes, border * sizeof(T));

      // Construct
      {
        char* ptr = buffer;
        unsigned nlines = 1;
        for (unsigned i = 0; i < dim-1; ++i)
          nlines *= shape[i];
        unsigned nelements = shape[dim-1];

        if (dim > 1) {
          for (unsigned i = 0; i < nlines; ++i, ptr += strides[dim-2]) {
            T* p_ = (T*) ptr;
            for (unsigned j = 0; j < nelements; ++j, ++p_)
              new (p_) T(v);
          }
        } else {
          T* p_ = (T*) ptr;
          for (unsigned j = 0; j < nelements; ++j, ++p_)
            new (p_) T(v);
        }
      }
    }


    template <typename T, unsigned dim>
    ndimage_data<T, dim>::~ndimage_data()
    {
      char* ptr = buffer;
      unsigned nlines = 1;
      for (unsigned i = 0; i < dim-1; ++i)
        nlines *= shape[i];
      unsigned nelements = shape[dim-1];
      if (dim == 1) {
        for (unsigned k = 0; k < nelements; ++k)
          ((T*)ptr + k)->~T();
        //    std::destroy( (T*)ptr, (T*) ptr + nelements);
      } else {
        for (unsigned i = 0; i < nlines; ++i, ptr += strides[dim-2])
          for (unsigned k = 0; k < nelements; ++k)
            ((T*)ptr + k)->~T();
        //std::destroy((T*)ptr, (T*)ptr + nelements);
      }
      mln::aligned_free(buffer);
    }
  }


  template <typename T, unsigned dim, typename E>
  ndimage_base<T,dim,E>::ndimage_base(unsigned border)
    : domain_ (), border_ (border), ptr_ (NULL)
  {
    for (unsigned i = 0; i < dim; ++i){
      mln_postcondition(domain_.pmin[i] == 0);
      mln_postcondition(domain_.pmax[i] == 0);
    }
  }

  template <typename T, unsigned dim, typename E>
  ndimage_base<T,dim,E>::ndimage_base(const domain_type& domain, unsigned border, T v)
    : domain_ (domain),
      border_ (border)
  {
    resize_(domain_, border, v);
  }

  template <typename T, unsigned dim, typename E>
  template <typename U, typename E2>
  ndimage_base<T,dim,E>::ndimage_base(const ndimage_base<U, dim, E2>& g, unsigned border, T v)
    : domain_ (g.domain()),
      border_ (border)
  {
    resize_(domain_, border_, v);
  }

  template <typename T, unsigned dim, typename E>
  template <typename U, typename E2>
  ndimage_base<T,dim,E>::ndimage_base(const ndimage_base<U, dim, E2>& g, mln::init)
    : domain_ (g.domain()),
      border_ (g.border())
  {
    resize_(domain_, border_, T());
  }

  template <typename T, unsigned dim, typename E>
  E
  ndimage_base<T,dim,E>::from_buffer_extern_(void* buffer,
                                             const domain_type& domain,
                                             const size_t* strides)
  {
    E image_;
    image_.domain_ = domain;
    image_.border_ = 0;

    MLN_EVAL_IF_DEBUG(image_.vbox_ = domain);
    std::copy(strides, strides + dim, image_.strides_.begin());

    if (image_.strides_[dim-1] % sizeof(T) != 0) {
      throw std::runtime_error("The padding of the image is not compatible with the size of the element.");
    }

    point<size_t, dim> sz = domain.shape();
    image_.m_ptr_origin = (T*)buffer;
    image_.m_index_strides[dim-1] = 1;
    image_.m_index_first = 0;
    image_.m_index_last = sz[dim-1] - 1;
    image_.ptr_ = (char*)buffer;
    image_.last_ = (char*)buffer + (sz[dim-1] - 1) * image_.strides_[dim-1];

    for (int i = dim-2; i >= 0; --i)
      {
        if (image_.strides_[i] % sizeof(T) != 0) {
          throw std::runtime_error("The padding of the image is not compatible with the size of the element.");
        }

        image_.m_index_strides[i] = image_.strides_[i] / sizeof(T);
        image_.m_index_last  += (sz[i] - 1) * image_.m_index_strides[i];
        image_.last_ += (sz[i] - 1) * image_.strides_[i];
      }

    return image_;
  }

  template <typename T, unsigned dim, typename E>
  E
  ndimage_base<T,dim,E>::from_buffer_copy_(void* buffer,
                                           const domain_type& domain,
                                           const size_t* strides)
  {
    E in = from_buffer_extern_(buffer, domain, strides);

    E out(domain);

    mln_pixter(pin, pout, in, out);
    mln_forall(pin, pout)
      pout->val() = pin->val();

    return out;
  }

  template <typename T, unsigned dim, typename E>
  E
  ndimage_base<T,dim,E>::from_buffer(void* buffer,
                                     const domain_type& domain,
                                     const size_t* strides,
                                     bool copy)
  {
    if (copy)
      return from_buffer_copy_(buffer, domain, strides);
    else
      return from_buffer_extern_(buffer, domain, strides);
  }



  template <typename T, unsigned dim, typename E>
  E
  ndimage_base<T,dim,E>::from_buffer(void* buffer,
                                     const domain_type& domain,
                                     bool copy)
  {
    point<size_t, dim> sz = domain.shape();
    size_t strides[dim];

    strides[dim-1] = sizeof(T);
    for (int i = dim-2; i >= 0; --i)
      strides[i] = strides[i+1] * sz[i+1];

    return from_buffer(buffer, domain, strides, copy);
  }

  template <typename T, unsigned dim, typename E>
  void
  ndimage_base<T,dim,E>::resize_(const domain_type& domain, unsigned border, T v)
  {
    site_type shp = domain.shape();
    point<size_t, dim> sz;
    MLN_EVAL_IF_DEBUG(vbox_ = domain);
    MLN_EVAL_IF_DEBUG(vbox_.pmin -= border);
    MLN_EVAL_IF_DEBUG(vbox_.pmax += border);
    sz = shp;

    // Compute strides size (in bytes)
    // The row stride is 16 bytes aligned
    data_.reset(new internal::ndimage_data<T, dim>(&(sz[0]), border, v));
    std::copy(data_->strides, data_->strides + dim, strides_.begin());

    // Compute pointer at (0,0)
    m_ptr_origin = (T*)data_->buffer;
    m_index_strides[dim-1] = 1;
    m_index_first = border_;
    m_index_last  = border_ + sz[dim-1] - 1;
    ptr_  = data_->buffer + border * strides_[dim-1];
    last_ = data_->buffer + (border + sz[dim-1] - 1) * strides_[dim-1];
    for (int i = dim-2; i >= 0; --i)
      {
        m_index_strides[i] = strides_[i] / sizeof(T);
        m_index_first += border_ * m_index_strides[i];
        m_index_last  += (border_ + sz[i] - 1) * m_index_strides[i];
        ptr_ += border * strides_[i];
        last_ += (border + sz[i] - 1) * strides_[i];
      }
  }


  template <typename T, unsigned dim, typename E>
  void
  ndimage_base<T,dim,E>::resize(const domain_type& domain, unsigned border, T v)
  {
    if (not (domain == domain_) or (int)border != border_)
      {
        domain_ = domain;
        border_ = border;
        resize_(domain_, border, v);
      }
    else
      {
        std::fill((T*) ptr_, ((T*)last_) + 1, v);
      }
  }

  template <typename T, unsigned dim, typename E>
  inline
  T&
  ndimage_base<T,dim,E>::at (const site_type& p)
  {
    mln_precondition(vbox_.has(p));

    char* ptr = ptr_;
    for (unsigned i = 0; i < dim; ++i)
      ptr += (p[i] - domain_.pmin[i]) * strides_[i];
    return * (T*)ptr;
  }


  template <typename T, unsigned dim, typename E>
  inline
  const T&
  ndimage_base<T,dim,E>::at (const site_type& p) const
  {
    mln_precondition(vbox_.has(p));

    char* ptr = ptr_;
    for (unsigned i = 0; i < dim; ++i)
      ptr += (p[i] - domain_.pmin[i]) * strides_[i];
    return * (const T*)ptr;
  }


  template <typename T, unsigned dim, typename E>
  inline
  T&
  ndimage_base<T,dim,E>::operator() (const site_type& p)
  {
    mln_precondition(domain_.has(p));
    return at(p);
  }

  template <typename T, unsigned dim, typename E>
  inline
  const T&
  ndimage_base<T,dim,E>::operator() (const site_type& p) const
  {
    mln_precondition(domain_.has(p));
    return at(p);
  }

  template <typename T, unsigned dim, typename E>
  inline
  T&
  ndimage_base<T,dim,E>::operator[] (size_type i)
  {
    mln_precondition(vbox_.has(point_at_index(i)));
    return *(m_ptr_origin + i);
  }

  template <typename T, unsigned dim, typename E>
  inline
  const T&
  ndimage_base<T,dim,E>::operator[] (size_type i) const
  {
    mln_precondition(vbox_.has(point_at_index(i)));
    return *(m_ptr_origin + i);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::pixel_type
  ndimage_base<T,dim,E>::pixel_at(const point_type& p)
  {
    pixel_type pix;
    pix.point_ = p;
    pix.ima_  = (E*)this;
    pix.ptr_ = (char*) & (operator () (p));
    return pix;
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::const_pixel_type
  ndimage_base<T,dim,E>::pixel_at(const point_type& p) const
  {
    const_pixel_type pix((const E*) this);
    pix.point_ = p;
    pix.ptr_ = (char*) & (operator () (p));
    return pix;
  }


  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::pixel_type
  ndimage_base<T,dim,E>::pixel(const point_type& p)
  {
    mln_precondition(domain_.has(p));
    return pixel_at(p);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::const_pixel_type
  ndimage_base<T,dim,E>::pixel(const point_type& p) const
  {
    mln_precondition(domain_.has(p));
    return pixel_at(p);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::size_type
  ndimage_base<T,dim,E>::index_of_point(const point_type& p) const
  {
    std::size_t idx = m_index_first;
    point_type  q = p - domain_.pmin;
    for (unsigned i = 0; i < dim; ++i)
      idx += q[i] * m_index_strides[i];
    return idx;
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::point_type
  ndimage_base<T,dim,E>::point_at_index(size_type idx) const
  {
    int k = idx;
    int kpmin = m_index_first;
    point_type p = point_type ();

    for (unsigned i = 0; i < dim; ++i) {
      std::div_t off = std::div((int)kpmin,  (int)m_index_strides[i]);
      std::div_t res = std::div((int)k,  (int)m_index_strides[i]);
      p[i] = res.quot - off.quot + domain_.pmin[i];
      k = res.rem;
      kpmin = off.rem;
    }
    mln_postcondition(vbox_.has(p));
    return p;
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::difference_type
  ndimage_base<T,dim,E>::delta_index(const point_type& p) const
  {
    difference_type idx = 0;
    for (unsigned i = 0; i < dim; ++i)
      idx += p[i] * m_index_strides[i];
    return idx;
  }

  // template <typename T, unsigned dim, typename E>
  // inline
  // T&
  // ndimage_base<T,dim,E>::element (difference_type n)
  // {
  //   return *reinterpret_cast<T*>(ptr_+n);
  // }

  // template <typename T, unsigned dim, typename E>
  // inline
  // const T&
  // ndimage_base<T,dim,E>::element (difference_type n) const
  // {
  //   return *reinterpret_cast<const T*>(ptr_+n);
  // }


  template <typename T, unsigned dim, typename E>
  const size_t*
  ndimage_base<T,dim,E>::strides () const
  {
    return &strides_[0];
  }


  // template <typename T, unsigned dim, typename E>
  // ptrdiff_t
  // ndimage_base<T,dim,E>::offset (point_type dp) const
  // {
  //   ptrdiff_t x = 0;
  //   for (unsigned i = 0; i < dim; ++i)
  //     x += p[i] * strides_[i];
  //   return x;
  // }

  template <typename T, unsigned dim, typename E>
  inline
  const typename ndimage_base<T,dim,E>::domain_type&
  ndimage_base<T,dim,E>::domain () const
  {
    return domain_;
  }

  /* -- Value range -- */

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::const_value_range
  ndimage_base<T,dim,E>::values () const
  {
    return const_value_range(exact(*this));
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::value_range
  ndimage_base<T,dim,E>::values ()
  {
    return value_range(*this);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::const_pixel_range
  ndimage_base<T,dim,E>::pixels () const
  {
    return const_pixel_range(*this);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::pixel_range
  ndimage_base<T,dim,E>::pixels ()
  {
    return pixel_range(*this);
  }

  template <typename T, unsigned dim, typename E>
  inline
  typename ndimage_base<T,dim,E>::extension_type
  ndimage_base<T,dim,E>::extension () const
  {
    return extension_type(ptr_, &strides_[0], domain_.shape(), border_);
  }

  template <unsigned d, typename T1, typename E1, typename T2, typename E2>
  inline
  bool are_indexes_compatible(const ndimage_base<T1, d, E1>& self,
                              const ndimage_base<T2, d, E2>& other)
  {
    return
      (self.index_of_point(self.domain().pmin) ==
       other.index_of_point(other.domain().pmin)) and
      (self.m_index_strides == other.m_index_strides);
  }



  /******************************************/
  /****     Extension implementation     ****/
  /******************************************/


  namespace internal
  {

    template <typename T, unsigned dim>
    ndimage_extension<T, dim>::ndimage_extension(char* ptr, const std::size_t* strides, const P& shp, int border)
      : m_shp (shp), m_border (border)
    {
      m_ptr = ptr;
      for (unsigned i = 0; i < dim; ++i) {
        m_strides[i] = strides[i];
        m_ptr -= m_strides[i] * border;
      }
    }

    template <typename T, unsigned dim>
    void
    ndimage_extension<T, dim>::fill(const T& v)
    {
      _fill<0>(m_ptr, v);
    }

    template <typename T, unsigned dim>
    void
    ndimage_extension<T, dim>::mirror()
    {
      _mirror<0>(m_ptr);
    }

    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d < dim)>::type
    ndimage_extension<T, dim>::_mirror(char* ptr)
    {
      char* pori = ptr + m_border * m_strides[d];
      char* p = pori;
      for (int i = 0; i < m_shp[d]; ++i, p += m_strides[d])
        _mirror<d+1>(p);

      char* src1 = pori;
      char* dest1 = pori - m_strides[d];
      char* dest2 = pori + m_shp[d] * m_strides[d];
      char* src2 = dest2 - m_strides[d];

      for (int i = 0; i < m_border; ++i)
        {
          _copy_line<d+1>(src1, dest1);
          _copy_line<d+1>(src2, dest2);
          src1 += m_strides[d];
          dest1 -= m_strides[d];
          src2 -= m_strides[d];
          dest2 += m_strides[d];
        }
    }

    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d == dim)>::type
    ndimage_extension<T, dim>::_mirror(char* ptr)
    {
      (void) ptr;
    }


    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d < dim)>::type
    ndimage_extension<T, dim>::_copy_line(char* src, char* dest)
    {
      for (int k = 0; k < m_shp[d] + 2 * m_border; ++k)
        {
          _copy_line<d+1>(src, dest);
          src += m_strides[d];
          dest += m_strides[d];
        }
    }

    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d == dim)>::type
    ndimage_extension<T, dim>::_copy_line(char* src, char* dest)
    {
      *(reinterpret_cast<T*>(dest)) = *(reinterpret_cast<T*>(src));
    }


    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d < dim)>::type
    ndimage_extension<T, dim>::_fill(char* ptr, const T& v)
    {
      for (int i = 0; i < m_border; ++i, ptr += m_strides[d])
        _fillall<d+1>(ptr, v);

      for (int i = 0; i < m_shp[d]; ++i, ptr += m_strides[d])
        _fill<d+1>(ptr, v);

      for (int i = 0; i < m_border; ++i, ptr += m_strides[d])
        _fillall<d+1>(ptr, v);
    }

      template <typename T, unsigned dim>
      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      ndimage_extension<T, dim>::_fill(char* ptr, const T& v)
      {
        (void) ptr;
        (void) v;
      }

    template <typename T, unsigned dim>
    template <unsigned d>
    typename std::enable_if<(d < dim)>::type
    ndimage_extension<T, dim>::_fillall(char* ptr, const T& v)
    {
      for (int i = 0; i < m_shp[d] + (2*m_border); ++i, ptr += m_strides[d])
        _fillall<d+1>(ptr, v);
    }

      template <typename T, unsigned dim>
      template <unsigned d>
      typename std::enable_if<(d == dim)>::type
      ndimage_extension<T, dim>::_fillall(char* ptr, const T& v)
      {
        *(reinterpret_cast<T*>(ptr)) = v;
      }

  }

} // end of namespace mln



#endif // !MLN_CORE_IMAGE_NDIMAGE_HH
