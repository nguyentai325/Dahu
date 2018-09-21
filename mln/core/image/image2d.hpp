#ifndef MLN_CORE_IMAGE_IMAGE2D_HPP
# define MLN_CORE_IMAGE_IMAGE2D_HPP

# include <mln/core/image/ndimage.hpp>
# include <mln/core/ch_value.hpp>

namespace mln {

  // FWD declaration
  template <typename T> struct image2d;


  // Specialization of traits
  template <typename T>
  struct image_traits< image2d<T> > : image_traits< ndimage_base<T, 2, image2d<T> > >
  {
  };

  template <typename T, typename V>
  struct image_ch_value< image2d<T>, V >
  {
    typedef image2d<V> type;
  };

  ///
  /// \brief The standard 2D image type.
  ///
  /// \p image2d is a type to represent 2d image. They are implemented with contiguous
  /// data buffer to allow pointer arithmetic. However data are not contiguous in this buffer
  /// since each row is aligned on a 128-bits boundary to get performance. \p image2d are soft image
  /// since two images can share the same buffer of data to avoid uncessary copies.
  ///
  /// \p image2d models a Writable Point-Accessible Image concept.
  /// See \see image2d \see image3d
  template <typename T>
  struct image2d : ndimage_base<T, 2, image2d<T> >
  {
  private:
    typedef ndimage_base<T, 2, image2d<T> > base_type;
    typedef ndimage_base<T, 2, image2d<T> > base;

  public:
    typedef typename base_type::domain_type domain_type;


    explicit image2d (unsigned border = 3)
      : ndimage_base<T,2,image2d<T> > (border)
    {
    }

    explicit image2d(const domain_type& domain, unsigned border = 3, const T& init = T())
      : ndimage_base<T,2, image2d<T> >(domain, border, init)
    {
    }

    image2d(short nrows, short ncols, unsigned border = 3)
      : ndimage_base<T,2, image2d<T> >( (box<short,2>){{0,0},{nrows, ncols}}, border)
    {
    }

    /// \brief Initialization constructor
    /// \{
    template <typename U>
    image2d(const image2d<U>& f, mln::init)
      : base(f, mln::init() )
    {
    }

    template <typename U>
    image2d(const image2d<U>& f, unsigned border, const T& v = T())
      : base(f, border, v)
    {
    }
    /// \}



    image2d(std::initializer_list< std::initializer_list<T> > l, unsigned border = 3)
    {
      mln_precondition(l.size() != 0);
      mln_precondition((l.begin())->size() != 0);

      short nr = l.size();
      short nc = (l.begin()->size());
      this->resize(domain_type{{0,0}, {nr,nc}}, border);

      mln_iter(v, this->values());
      v.init();
      for (const std::initializer_list<T>& row : l)
        {
          mln_assertion(row.size() == (unsigned)nc);
          for (const T* val = row.begin(); val != row.end(); ++val, v.next())
            *v = *val;
        }
    }

    using base::at;

    T& at(short nrows, short ncols)
    {
      return base::at(typename base::point_type(nrows, ncols));
    }

    const T& at(short nrows, short ncols) const
    {
      return base::at(typename base::point_type(nrows, ncols));
    }

    unsigned nrows() const { return this->domain_.pmax[0] - this->domain_.pmin[0]; }
    unsigned ncols() const { return this->domain_.pmax[1] - this->domain_.pmin[1]; }

  };



} // end of namespace mln



#endif //!MLN_CORE_IMAGE_IMAGE2D_HPP
