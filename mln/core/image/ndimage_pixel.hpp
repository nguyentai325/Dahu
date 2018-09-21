#ifndef MLN_CORE_IMAGE_NDIMAGE_PIXEL_HPP
# define MLN_CORE_IMAGE_NDIMAGE_PIXEL_HPP

# include <mln/core/image/ndimage.hpp>

namespace mln
{

  template <typename T, unsigned dim, typename I>
  struct ndimage_pixel : Pixel< ndimage_pixel<T, dim, I> >
  {
    typedef mln::point<short, dim>              site_type;
    typedef site_type				point_type;
    typedef typename std::remove_const<T>::type value_type;
    typedef T&                                  reference;
    typedef ptrdiff_t                           distance_type;
    typedef I                                   image_type;
    typedef typename I::size_type		size_type;

    ndimage_pixel() : ima_ (NULL) {}

    ndimage_pixel(image_type* ima) : ima_ (ima) {}


    /// \brief Copy constructor.
    ndimage_pixel(const ndimage_pixel& pix) = default;

    /// \brief Interopability constructor from the pixel type to the const one.
    template <typename U, typename Other>
    ndimage_pixel(const ndimage_pixel<U, dim, Other>& pix,
                  typename std::enable_if< std::is_convertible<U*,T*>::value, void*>::type = NULL);


    /// \brief Assignment operator
    ndimage_pixel&      operator= (const ndimage_pixel& other);


    reference	        val() const;
    site_type           point()  const;
    site_type           site()  const;
    distance_type	offset() const;
    image_type&         image() const;
    size_type		index() const { return index_; }

    template <typename, unsigned, typename> friend struct ndimage_pixel;
    template <typename, unsigned, typename> friend struct ndimage_base;

  private:
    friend struct internal::iterator_core_access;

    char* &		get_value() { return ptr_; }
    site_type&		get_point() { return point_; }
    size_t&		get_index() { return index_; }
    const char*         get_value() const { return ptr_; }
    const site_type&    get_point() const { return point_; }
    size_t		get_index() const { return index_; }

    template <typename U, typename Other>
    typename std::enable_if< std::is_convertible<U*,T*>::value, bool>::type
    equal(const ndimage_pixel<U, dim, Other>& other) const
    {
      return this->ptr_ == other.ptr_;
    }


  private:
    image_type*         ima_;
    site_type           point_;
    char*               ptr_;
    size_t		index_;
  };

/******************************************/
/****          Implementation          ****/
/******************************************/


  template <typename T, unsigned dim, typename I>
  template <typename U, typename Other>
  inline
  ndimage_pixel<T, dim, I>::ndimage_pixel(const ndimage_pixel<U, dim, Other>& pix,
                                          typename std::enable_if< std::is_convertible<U*,T*>::value, void*>::type)
    : ima_ (pix.ima_),
      point_ (pix.point_),
      ptr_ (pix.ptr_),
      index_ (pix.index_)
  {
  }

  template <typename T, unsigned dim, typename I>
  inline
  ndimage_pixel<T, dim, I>&
  ndimage_pixel<T, dim, I>::operator= (const ndimage_pixel& other)
  {
    //mln_precondition(other.ima_ == ima_);
    ima_ = other.ima_;
    point_ = other.point_;
    ptr_ = other.ptr_;
    index_ = other.index_;
    return *this;
  }

  template <typename T, unsigned dim, typename I>
  inline
  typename ndimage_pixel<T,dim,I>::reference
  ndimage_pixel<T,dim,I>::val() const
  {
    return *reinterpret_cast<T*>(ptr_);
  }

  template <typename T, unsigned dim, typename I>
  inline
  typename ndimage_pixel<T,dim,I>::site_type
  ndimage_pixel<T,dim,I>::point() const
  {
    return point_;
  }

  template <typename T, unsigned dim, typename I>
  inline
  typename ndimage_pixel<T,dim,I>::site_type
  ndimage_pixel<T,dim,I>::site() const
  {
    return point_;
  }

  template <typename T, unsigned dim, typename I>
  inline
  typename ndimage_pixel<T,dim,I>::distance_type
  ndimage_pixel<T,dim,I>::offset() const
  {
    mln_precondition(ima_ != NULL && ptr_ != NULL);
    return static_cast<char*>(ptr_) - ima_->ptr_;
  }

  template <typename T, unsigned dim, typename I>
  inline
  typename ndimage_pixel<T,dim,I>::image_type&
  ndimage_pixel<T,dim,I>::image() const
  {
    mln_precondition(ima_ != NULL);
    return *ima_;
  }


}

#endif
