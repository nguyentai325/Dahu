#ifndef MORPHER_BASE_HPP
# define MORPHER_BASE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/morphers/details/morpher_core_access.hpp>
# include <mln/core/pixel_utility.hpp>

//
// This interface aims at easing the creation of image morpher.
// To create a simple morpher just inherits from morpher_base,
// and define the way the interface can access the underlying
// morphed image, either by:
// * overriding `image_type& get_morphed() const` and `const image_type& get_morphed() const` methods
// * declaring a variable `image_type m_ima`
//
// The previous accessors shall be made private and morpher_core_access
// has to be a friend class to allow the base interface implementation to
// access them.
//
//
//
// The following shows an example of a perfect forwarding image morpher.
// It does nothing but forwarding the method calls.
//
// template <class I>
// struct mymorpher : morpher_base<mymorpher, I>
// {
//   friend struct morpher_core_access;
//
//   mymorpher(I& image) : m_ima(image) {}
//   private:
//      I& m_ima;
// };
//
// and inherits the traits as well
//
// template <class I>
// struct image_traits<mymorpher<I>> : morpher_base_traits<I>
// {
// };
//
// Doing so results in a type that totally has the same
// methods than I (methods that are part of the Image interface)
// and with the same traits. Now you can override the methods
// to have a morpher that actually do something.


namespace mln
{

  /// \brief Base class for morpher that acts like a mixin w.r.t to I
  ///
  ///
  template <typename Derived, typename I, typename P = mln_point(I), typename V = mln_value(I)>
  struct morpher_base;

  /// \brief Base class for traits
  template <class I>
  struct morpher_base_traits : image_traits<I>
  {
  };

  /******************************************************/
  /**                HELPER MACROS                     **/
  /******************************************************/

  // Consider perfect forwarding here.

# define MLN_IMMORPHER_FORWARD_IF_0_(COND, F, RETURN, CV)       \
  typename std::enable_if<COND, RETURN>::type                   \
  F() CV                                                        \
  {                                                             \
    return morpher_core_access::get_ima(this).F();              \
  }

# define MLN_IMMORPHER_FORWARD_IF_1_(COND, F, RETURN, ARG, CV)   \
  typename std::enable_if<COND, RETURN>::type                   \
  F(ARG x) CV                                                   \
  {                                                             \
    return morpher_core_access::get_ima(this).F(x);             \
  }

# define MLN_IMMORPHER_FORWARD_0(F, RETURN)          MLN_IMMORPHER_FORWARD_IF_0_(true, F, RETURN, )
# define MLN_IMMORPHER_FORWARD_CONST_0(F, RETURN)    MLN_IMMORPHER_FORWARD_IF_0_(true, F, RETURN, const)
# define MLN_IMMORPHER_FORWARD_IF_0(COND, F, RETURN) MLN_IMMORPHER_FORWARD_IF_0_(COND, F, RETURN, )
# define MLN_IMMORPHER_FORWARD_IF_CONST_0(F, RETURN) MLN_IMMORPHER_FORWARD_IF_0_(COND, F, RETURN, const)

# define MLN_IMMORPHER_FORWARD_1(F, RETURN, ARG)                MLN_IMMORPHER_FORWARD_IF_1_(true, F, RETURN, ARG,)
# define MLN_IMMORPHER_FORWARD_CONST_1(F, RETURN, ARG)          MLN_IMMORPHER_FORWARD_IF_1_(true, F, RETURN, ARG, const)
# define MLN_IMMORPHER_FORWARD_IF_1(COND, F, RETURN, ARG)       MLN_IMMORPHER_FORWARD_IF_1_(COND, F, RETURN, ARG,)
# define MLN_IMMORPHER_FORWARD_IF_CONST_1(COND, F, RETURN, ARG) MLN_IMMORPHER_FORWARD_IF_1_(COND, F, RETURN, ARG, const)

  /***********************/
  /**  Implementation   **/
  /***********************/



  //
  // We define in this namespace some Mixin that will delegate the method calls
  // to the morphed image w.r.t. to the image categories.
  //
  namespace impl
  {
    /// \brief Default implementation for morphing an accessible image
    template <typename Derived, typename I, bool _is_accessible = image_traits<Derived>::accessible::value>
    struct morpher_accessible
    {
      friend struct mln::morpher_core_access;
    };

    template <typename Derived, typename I>
    struct morpher_accessible<Derived, I, true>
    {
    private:
      friend  struct mln::morpher_core_access;
      typedef I         image_t;
      typedef Derived   derived_t;

    public:
      typedef mln_point(I)                              point_type;
      typedef mln_point(I)                              site_type;
      typedef typename image_reference<I>::type		reference;
      typedef typename image_const_reference<I>::type	const_reference;

      MLN_IMMORPHER_FORWARD_1(operator(), reference, const point_type&);
      MLN_IMMORPHER_FORWARD_CONST_1(operator(), const_reference, const point_type&);
      MLN_IMMORPHER_FORWARD_1(at, reference, const point_type&);
      MLN_IMMORPHER_FORWARD_CONST_1(at, const_reference, const point_type&);
    };

    /// \brief Default implementation for morphing an accessible image
    template <typename Derived, typename I, bool _is_indexable = image_traits<Derived>::indexable::value>
    struct morpher_indexable
    {
    };

    template <typename Derived, typename I>
    struct morpher_indexable<Derived, I, true>
    {
    private:
      friend  struct mln::morpher_core_access;
      typedef I						image_t;
      typedef Derived                                   derived_t;
      typedef typename I::point_type			point_type;

    public:
      typedef typename I::size_type			size_type;
      typedef typename I::difference_type		difference_type;
      typedef typename image_reference<I>::type		reference;
      typedef typename image_const_reference<I>::type	const_reference;

      MLN_IMMORPHER_FORWARD_1(operator[], reference, size_type);
      MLN_IMMORPHER_FORWARD_CONST_1(operator[], const_reference, size_type);
      MLN_IMMORPHER_FORWARD_CONST_1(index_of_point, size_type, const point_type&);
      MLN_IMMORPHER_FORWARD_CONST_1(delta_index, difference_type, const point_type&);
      MLN_IMMORPHER_FORWARD_CONST_1(point_at_index, point_type, size_type);
    };
  }


  template <typename Derived, typename I, typename P, typename V>
  struct morpher_base : image_base<Derived, P, V>,
    impl::morpher_accessible<Derived, typename std::remove_reference<I>::type >,
    impl::morpher_indexable<Derived, typename std::remove_reference<I>::type >
  {
  private:
    friend  struct mln::morpher_core_access;
    typedef typename std::remove_reference<I>::type image_t;
    typedef Derived                                 derived_t;

  public:
    typedef typename image_t::value_type	value_type;
    typedef typename image_t::reference         reference;
    typedef typename image_t::const_reference   const_reference;
    typedef typename image_t::site_type		site_type;
    typedef typename image_t::point_type	point_type;



    typedef typename image_t::domain_type                       domain_type;
    typedef typename image_value_range<image_t>::type		value_range;
    typedef typename image_const_value_range<image_t>::type	const_value_range;
    typedef typename image_pixel_range<image_t>::type		pixel_range;
    typedef typename image_const_pixel_range<image_t>::type	const_pixel_range;


    auto domain() const
      -> decltype(morpher_core_access::get_ima(this).domain())
    {
      return morpher_core_access::get_ima(this).domain();
    }


    MLN_IMMORPHER_FORWARD_0(values, value_range);
    MLN_IMMORPHER_FORWARD_0(pixels, pixel_range);
    MLN_IMMORPHER_FORWARD_CONST_0(values, const_value_range);
    MLN_IMMORPHER_FORWARD_CONST_0(pixels, const_pixel_range);

  protected:
    image_t&		get_morphed() { return morpher_core_access::get_ima_(this); }
    const image_t&	get_morphed() const { return morpher_core_access::get_ima_(this); }
  };


  // template <class Pixel, class Morpher>
  // struct morpher_default_pixel
  //   : morpher_pixel_base< morpher_default_pixel<Pixel, Morpher>,
  //                         Pixel, Morpher >
  // {
  //   morpher_default_pixel() = default;
  //   morpher_default_pixel(const morpher_default_pixel&) = default;

  //   morpher_default_pixel(const Pixel& pix, Mopher& morpher)
  //     : morpher_pixel_base (morpher),
  //       m_pix (pix)
  //   {
  //   }

  //   template <class Pixel2, class Morpher2>
  //   morpher_default_pixel(const morpher_default_pixel<Pixel2, Morpher2>& other,
  //                         typename std::enable_if<std::is_convertible<Morpher2*, Morpher*>::value and
  //                                                 std::is_convertible<Pixel2, Pixel>::value>::type* = NULL)
  //   : morpher_pixel_base(other.image()),
  //     m_pix(other.m_pix)
  //   {
  //   }

  // private:
  //   Pixel m_pix;
  // };

}

# undef MLN_IMMORPHER_FORWARD_0_
# undef MLN_IMMORPHER_FORWARD_0
# undef MLN_IMMORPHER_FORWARD_CONST_0
# undef MLN_IMMORPHER_FORWARD_IF_0
# undef MLN_IMMORPHER_FORWARD_IF_CONST_0
# undef MLN_IMMORPHER_FORWARD_1_
# undef MLN_IMMORPHER_FORWARD_1
# undef MLN_IMMORPHER_FORWARD_CONST_1
# undef MLN_IMMORPHER_FORWARD_IF_1
# undef MLN_IMMORPHER_FORWARD_IF_CONST_1


#endif // ! MORPHER_BASE_HPP
