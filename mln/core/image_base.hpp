#ifndef MLN_CORE_IMAGE_IMAGE_HPP
# warning "Do not include this file directly but <mln/core/image/image.hpp>"
# include <mln/core/image/image.hpp>
#endif

#ifndef MLN_CORE_IMAGE_BASE_HPP
# define MLN_CORE_IMAGE_BASE_HPP

/// \file

# include <mln/core/image/image.hpp>
# include <cstddef>

namespace mln
{

  struct undef {};

  template <typename Derived,
            typename Point,
            typename Value,
            typename Pixel = undef,
            typename ConstPixel = undef,
            typename Reference = Value&,
            typename ConstReference = const Value&,
            typename Pointer = Value*,
            typename ConstPointer = const Value*,
            typename SizeType = std::size_t,
            typename Difference = std::ptrdiff_t
            >
  struct image_base : mln::Image<Derived>
  {
    typedef Point point_type;
    typedef Point site_type;
    typedef Pixel pixel_type;
    typedef ConstPixel const_pixel_type;

    typedef Value value_type;
    typedef Reference reference;
    typedef ConstReference const_reference;
    typedef Pointer pointer;
    typedef ConstPointer const_pointer;

  };

  /******************************************/
  /****   Image Documentation Methods   *****/
  /******************************************/
#ifdef MLN_DOXYGEN

  /// \brief Documentation class.
  struct image
  {
    /// \brief Type of point/site of the image.
    typedef implementation_defined site_type;

    /// \brief Alias for site_type
    typedef implementation_defined point_type;

    /// \brief The value type of the image.
    typedef implementation_defined value_type;

    /// \brief The pixel type of the image.
    typedef implementation_defined pixel_type;

    /// \brief The non-mutable pixel type of the image.
    typedef implementation_defined pixel_type;

    /// \brief Reference type to a pixel value.
    typedef implementation_defined reference;

    /// \brief Const reference type to a pixel value.
    typedef implementation_defined const_reference;

    /// \brief Unsigned integral type for indexes.
    typedef implementation_defined size_type;

    /// \brief Signed integral type to represent difference between indexes.
    typedef implementation_defined difference_type;


    /**
    * \brief Access specified element by its index
    *
    * Accessing pixel values through their indexes allows a faster
    * processing compared to point/site access (see operator()(const
    * site_type& p) const ).
    *
    * \param i index of the element
    * \return Reference to the requested element.
    */
    const_reference operator[] (size_type i) const;

    /**
    * \brief Access specified element by its point/site
    *
    * Access a pixel values through their point/site. Note that the
    * point must belong to the domain of the image. To access the
    * extension, use at(const site_type&).
    *
    * \param p location of the element
    * \return Reference to the requested element
    * \pre \p p must belongs to the domain
    *
    */
    const_reference operator() (const point_type& p) const;

    /**
    * \brief Access specified element by its point/site without bound checking
    *
    * This accessor allows to access a value which is in the extension
    * of the image, contrary to the operator () that allows to access
    * only values within the image domain.
    *
    * \param p location of the element
    * \return Reference to the requested element
    * \pre \p p must belongs to the domain.
    */
    const_reference at() (const point_type& p) const;

    /// \brief The domain type of the image.
    typedef implementation_defined domain_type;

    /// \brief A range type to traverse image values.
    typedef implementation_defined value_range;

    /// \brief A range type to traverse image values.
    typedef implementation_defined const_value_range;

    /// \brief A range type to traverse image values.
    typedef implementation_defined pixel_range;

    /// \brief A range type to traverse image pixels.
    typedef implementation_defined const_pixel_range;


    /**
     * \brief Return a reference to the domain of the image.
     *
     * The domain is a sorted range so that a forward traversal
     * sees the site set increasing.
     */
    const domain_type& domain() const;

    /// \brief Return a const range over the values of the image.
    const_value_range values() const;

    /// \brief Return a const range over the pixels of the image.
    const_pixel_range pixels() const;

    /// \brief Return a range over the values of the image.
    value_range values();

    /// \brief Return a range over the pixels of the image.
    pixel_range pixels();


  /// \brief Return the index of the given point \p p.
  size_type index_of_point(const point_type& p) const;

  /// \brief Return the index difference given a delta point \p p.
  difference_type delta_index(const point_type& p) const;

  /// \brief Return the point of the given index \p i.
  point_type point_at_index(size_type i) const;

  /// \brief Return the pixel of the given point \p p.
  const_pixel_type pixel(const point_type& p) const;

  /// \brief Return the pixel of the given point \p p without bound checking.
  const_pixel_type pixel_at(const point_type& p) const;


  /// \brief Return the extension of the image.
  extension_type extension() const;

  };
#endif

} // end of namespace mln

#endif //!MLN_CORE_IMAGE_BASE_HPP
