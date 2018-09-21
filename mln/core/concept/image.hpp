#ifndef MLN_CORE_CONCEPT_IMAGE_HPP
# define MLN_CORE_CONCEPT_IMAGE_HPP

# include <mln/core/concept/object.hpp>
# include <mln/core/concept/check.hpp>
# include <mln/core/concept/pixel.hpp>
# include <mln/core/concept/extension.hpp>

# include <mln/core/image_traits.hpp>

# include <boost/concept_check.hpp>
# include <boost/static_assert.hpp>


namespace mln {

  template <typename I> struct Image;
  template <typename I> struct Image_;
  template <typename I> struct ConcreteImage;
  template <typename I> struct IterableImage;
  template <typename I> struct IndexableImage;
  template <typename I> struct AccessibleImage;
  template <typename I> struct ImageWithExtension;

  template <typename I>
  struct Image : Object_<I>
  {
    BOOST_CONCEPT_USAGE(Image)
    {
      BOOST_CONCEPT_ASSERT((Image_<I>));
    }

  };

  template <typename I>
  struct Image_
  {
    typedef image_traits<I> traits;

    typedef typename traits::accessible   accessible;
    typedef typename traits::category     category;
    typedef typename traits::concrete     concrete;
    typedef typename traits::indexable    indexable;
    typedef typename traits::extension    extension;

    typedef typename I::value_type         value;
    typedef typename I::reference          reference;
    typedef typename I::const_reference    const_reference;
    typedef typename I::pixel_type         pixel;
    typedef typename I::const_pixel_type   const_pixel;


    typedef typename I::site_type          site_type;
    typedef typename I::point_type         point_type;
    typedef typename I::domain_type        domain_type;

    BOOST_CONCEPT_ASSERT((Pixel<pixel>));
    BOOST_CONCEPT_ASSERT((Pixel<const_pixel>));

    BOOST_CONCEPT_USAGE(Image_)
    {
      check(std::is_base_of< Image<I>, I > ());
      // FIXME: to relax, an image can return an rvalue for the domain
      // const domain_type& (I::*method) () const = &I::domain;
      // (void) method;

      domain_type dom = ima.domain();
      (void) dom;

      check(std::is_convertible<typename pixel::value_type, value> ());
      check(std::is_same<typename pixel::reference, reference> ());
      check(std::is_convertible<typename const_pixel::value_type, value> ());
      check(std::is_same<typename const_pixel::reference, const_reference> ());
      check(std::is_convertible<pixel, const_pixel> ());

      {
        mln_concrete(I)       __ima2 = imconcretize(ima);
        mln_ch_value(I, int)  __ima3 = imchvalue<int>(ima);
      }

      MLN_CONCEPT_ASSERT_IF(image_traits<I>::concrete::value,
                            ConcreteImage<I>);

      MLN_CONCEPT_ASSERT_IF((std::is_convertible<category, forward_image_tag>::value),
                            IterableImage<I>);

      MLN_CONCEPT_ASSERT_IF(image_traits<I>::indexable::value,
                            IndexableImage<I>);

      MLN_CONCEPT_ASSERT_IF(image_traits<I>::accessible::value,
                            AccessibleImage<I>);

      MLN_CONCEPT_ASSERT_IF(image_has_extension<I>::value,
                            ImageWithExtension<I>);

    }

  private:
    I ima;
  };

  template <typename I>
  struct ImageWithExtension
  {
    BOOST_CONCEPT_USAGE(ImageWithExtension)
    {
      typedef typename I::extension_type extension;
      BOOST_CONCEPT_ASSERT((Extension<extension>));
    }
  };


  template <typename I>
  struct ConcreteImage : Image_<I>
  {
  public:
    typedef typename image_traits<I>::concrete concrete;

    static_assert(concrete::value, "Image must be concrete");

    template <class J, bool has_border = image_has_border<I>::value>
    struct check_border
    {
      BOOST_CONCEPT_USAGE(check_border)
      {
        mln_value(J) val;
        J f = *((J*)0);
        J g(f, mln::init ());
        J h(f, f.border());
        J i(f, f.border(), val);
      }
    };

    template <class J>
    struct check_border<J, false>
    {
      BOOST_CONCEPT_USAGE(check_border)
      {
        mln_value(J) val;
        J f = *((J*)0);
        J g(f, mln::init ());
        J h(f, val);
      }
    };

    template <class J>
    struct check_indexable
    {
      BOOST_CONCEPT_USAGE(check_indexable)
      {
        typedef typename I::size_type size_type;
        void (I::*ptr) (size_type) = &I::reindex;
        (void) ptr;
      }
    };


    BOOST_CONCEPT_USAGE(ConcreteImage)
    {
      BOOST_CONCEPT_ASSERT((check_border<I>));
      MLN_CONCEPT_ASSERT_IF(image_traits<I>::indexable::value,
                            check_indexable<I>);
    }

  };

  template <typename I>
  struct AccessibleImage : Image_<I>
  {
  public:
    typedef typename image_traits<I>::accessible accessible;

    static_assert(accessible::value, "Image must be accessible");

    BOOST_CONCEPT_USAGE(AccessibleImage)
    {
      typedef typename I::point_type        point_type;
      typedef typename I::reference         reference;
      typedef typename I::const_reference   const_reference;
      typedef typename I::pixel_type        pixel_type;
      typedef typename I::const_pixel_type  const_pixel_type;

      reference (I::*ptr) (const point_type&) = &I::operator();
      const_reference (I::*ptr2) (const point_type&) const = &I::operator();
      reference (I::*ptr3) (const point_type&) = &I::at;
      const_reference (I::*ptr4) (const point_type&) const = &I::at;
      pixel_type (I::*ptr5) (const point_type&) = &I::pixel;
      const_pixel_type (I::*ptr6) (const point_type&) const = &I::pixel;

      (void) ptr; (void) ptr2; (void) ptr3; (void) ptr4;
      (void) ptr5; (void) ptr6;
    }
  };


  template <typename I>
  struct IndexableImage : Image_<I>
  {
  public:
    typedef typename image_traits<I>::indexable indexable;

    static_assert(indexable::value, "Image must be indexable");

    BOOST_CONCEPT_USAGE(IndexableImage)
    {
      typedef typename I::size_type	    size_type;
      typedef typename I::difference_type   difference_type;
      typedef typename I::reference         reference;
      typedef typename I::const_reference   const_reference;
      typedef typename I::point_type        point_type;

      reference (I::*ptr) (size_type) = &I::operator[];
      const_reference (I::*ptr2) (size_type) const = &I::operator[];
      size_type (I::*ptr3) (const point_type&) const = &I::index_of_point;
      point_type (I::*ptr4) (size_type) const = &I::point_at_index;
      difference_type (I::*ptr5) (const point_type&) const = &I::delta_index;

      (void) ptr; (void) ptr2; (void) ptr3; (void) ptr4; (void) ptr5;
    }
  };


  template <typename I>
  struct IterableImage : Image_<I>
  {
  public:
    typedef typename image_traits<I>::category category;

    static_assert(std::is_convertible<category, forward_image_tag>::value,
                  "Image category must be iterable");

    BOOST_CONCEPT_USAGE(IterableImage)
    {
      typedef typename I::value_type         value_type;
      typedef typename I::reference         reference;
      typedef typename I::const_reference   const_reference;
      typedef typename I::pixel_type        pixel_type;
      typedef typename I::const_pixel_type  const_pixel_type;

      typedef typename I::value_range               value_range;
      typedef typename I::const_value_range         const_value_range;
      typedef typename I::pixel_range               pixel_range;
      typedef typename I::const_pixel_range         const_pixel_range;


      value_range       (I::*ptr1) () = &I::values;
      const_value_range (I::*ptr2) () const = &I::values;
      pixel_range       (I::*ptr3) () = &I::pixels;
      const_pixel_range (I::*ptr4) () const = &I::pixels;

      (void) ptr1;
      (void) ptr2;
      (void) ptr3;
      (void) ptr4;

      typedef typename value_range::iterator       value_iterator;
      typedef typename const_value_range::iterator const_value_iterator;
      typedef typename pixel_range::iterator       pixel_iterator;
      typedef typename const_pixel_range::iterator const_pixel_iterator;

       check(std::is_convertible<typename value_iterator::value_type, value_type> ());
      // "Iterator's value type is expected to be the image value type.");
       check(std::is_same<typename value_iterator::reference, reference> ());
      // "Iterator's reference is expected to be the image reference");
       check(std::is_convertible<typename const_value_iterator::value_type, value_type> ());
      // "Iterator's value type is expected to be the image value type");
      check(std::is_same<typename const_value_iterator::reference, const_reference> ());
      // "Iterator's reference type is expected to be the image const reference");

      //check(std::is_const<typename std::remove_reference<const_reference>::type > ());

      check(std::is_convertible<typename pixel_iterator::reference, pixel_type> ());
      check(std::is_same<typename pixel_iterator::value_type, pixel_type> ());
      // "Pixel Iterator's value type is expected to be the image pixel type");

      check(std::is_convertible<typename const_pixel_iterator::reference, const_pixel_type> ());
      check(std::is_same<typename const_pixel_iterator::value_type, const_pixel_type> ());
      // "Pixel Iterator's value type is expected to be the image const pixel type");

      check(std::is_same<typename pixel_type::image_type, I> ());
      check(std::is_same<typename const_pixel_type::image_type, const I> ());

      check(std::is_convertible<value_iterator, const_value_iterator> ());
      check(std::is_convertible<pixel_iterator, const_pixel_iterator> ());
    }
  };




} // end of namespace mln


#endif //!MLN_CORE_CONCEPT_IMAGE_HPP
