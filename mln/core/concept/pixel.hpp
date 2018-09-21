#ifndef MLN_CORE_CONCEPT_PIXEL_HPP
# define MLN_CORE_CONCEPT_PIXEL_HPP

# include <mln/core/concept/check.hpp>
# include <boost/concept_check.hpp>

namespace mln
{

  template <typename Pix>
  struct Pixel
  {

    BOOST_CONCEPT_USAGE(Pixel)
    {
      check(std::is_base_of<Pixel, Pix> ());

      typedef typename Pix::value_type        value_type;
      typedef typename Pix::reference         reference;
      typedef typename Pix::point_type        point_type;
      typedef typename Pix::site_type         site_type;
      typedef typename Pix::image_type        image_type;

      reference   (Pix::*method1) () const = &Pix::val;
      point_type  (Pix::*method2) () const = &Pix::point;
      site_type   (Pix::*method3) () const = &Pix::site;
      image_type& (Pix::*method4) () const = &Pix::image;

      (void) method1;
      (void) method2;
      (void) method3;
      (void) method4;

      // site() and point() aliase each others
      check(std::is_same<point_type, site_type> ());

      // value_type should not be a reference
      check_false(std::is_reference<value_type> ());
    }

  };

  template <class P, class V = int, class I = void*, class Ref = void>
  struct pixel_archetype : Pixel< pixel_archetype<P, V, I> >
  {
    struct default_reference {
      operator V() const { return std::declval<V>(); }
    };

    typedef P point_type;
    typedef P site_type;
    typedef V value_type;
    typedef I image_type;
    typedef typename std::conditional<std::is_same<Ref, void>::value,
				      default_reference,
				      Ref>::type reference;

    reference val() const { return make_object<reference>(); }
    P point() const { return make_object<P>(); }
    P site() const { return make_object<P>(); }
    I& image() const { return make_object<I&>(); }
  };


} // end of namespace mln

#endif //!MLN_CORE_CONCEPT_PIXEL_HPP
