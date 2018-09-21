#ifndef MLN_CORE_IMAGE_SUB_IMAGE_HPP
# define MLN_CORE_IMAGE_SUB_IMAGE_HPP

# include <type_traits>

# include <mln/core/image_base.hpp>
# include <mln/core/iterator/transform_iterator.hpp>
# include <mln/core/range/range_traits.hpp>
# include <mln/core/range/iterator_range.hpp>
# include <mln/core/range/filter.hpp>
# include <mln/core/pixel_utility.hpp>
# include <boost/optional.hpp>

namespace mln
{

  template <typename Image, typename Domain>
  struct sub_image;

  namespace internal
  {
    template <typename I>
    struct is_in_mask;
  }

  /// These overloads are worst-conversion possible
  /// \{

  // Ima | Domain
  template <typename I, typename Domain>
  sub_image<const I&, Domain>
  make_subimage(const Image<I>& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type* = NULL);

  template <typename I, typename Domain>
  sub_image<I&, Domain>
  make_subimage(Image<I>& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type* = NULL);

  template <typename I, typename Domain>
  sub_image<const I, Domain>
  make_subimage(const Image<I>&& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type* = NULL);

  template <typename I, typename Domain>
  sub_image<I, Domain>
  make_subimage(Image<I>&& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type* = NULL );


  // Ima | mask
  template <typename I, typename J>
  sub_image<const I&, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(const Image<I>& ima, const Image<J>& mask);

  template <typename I, typename J>
  sub_image<I&, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(Image<I>& ima, const Image<J>& mask);

  template <typename I, typename J>
  sub_image<const I, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(const Image<I>&& ima, const Image<J>& mask);

  template <typename I, typename J>
  sub_image<I, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(Image<I>&& ima, const Image<J>& mask);
  /// \}

  // Ima | {domain, mask}
  template <typename I, typename DomainOrMask>
  auto operator| (const Image<I>& ima, const DomainOrMask& other)
    -> decltype( make_subimage(exact(ima), other) );

  template <typename I, typename DomainOrMask>
  auto operator| (Image<I>& ima, const DomainOrMask& other)
    -> decltype( make_subimage(exact(ima), other) );

  template <typename I, typename DomainOrMask>
  auto operator| (const Image<I>&& ima, const DomainOrMask& other)
    -> decltype( make_subimage(move_exact(ima), other) );

  template <typename I, typename DomainOrMask>
  auto operator| (Image<I>&& ima, const DomainOrMask& other)
    -> decltype( make_subimage(move_exact(ima), other) );


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  /// These overloads are worst-conversion possible
  /// \{

  // Ima | Domain
  template <typename I, typename Domain>
  sub_image<const I&, Domain>
  make_subimage(const Image<I>& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type*)
  {
    return sub_image<const I&, Domain>(exact(ima), domain);
  }

  template <typename I, typename Domain>
  sub_image<I&, Domain>
  make_subimage(Image<I>& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type*)
  {
    return sub_image<I&, Domain>(exact(ima), domain);
  }

  template <typename I, typename Domain>
  sub_image<const I, Domain>
  make_subimage(const Image<I>&& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type*)
  {
    return sub_image<I, Domain>(move_exact(ima), domain);
  }

  template <typename I, typename Domain>
  sub_image<I, Domain>
  make_subimage(Image<I>&& ima, const Domain& domain, typename std::enable_if< not is_a<Domain, Image>::value >::type*)
  {
    return sub_image<I, Domain>(move_exact(ima), domain);
  }

  // Ima | Mask
  namespace internal
  {
    template <typename I>
    struct is_in_mask
    {
      typedef typename I::point_type P;

      is_in_mask(const I& ima)
      : m_mask(ima)
      {
      }

      bool operator () (const P& p) const
      {
	return m_mask(p);
      }

    private:
      const I& m_mask;
    };
  }


  template <typename I, typename J>
  sub_image<const I&, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(const Image<I>& ima, const Image<J>& mask)
  {
    static_assert(std::is_convertible<mln_value(J), bool>::value, "J's value type must be convertible to bool.");
    return make_subimage(exact(ima), rng::filter(exact(mask).domain(), internal::is_in_mask<J> (exact(mask))));
  }

  template <typename I, typename J>
  sub_image<I&, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(Image<I>& ima, const Image<J>& mask)
  {
    static_assert(std::is_convertible<mln_value(J), bool>::value, "J's value type must be convertible to bool.");
    return make_subimage(exact(ima), rng::filter(exact(mask).domain(), internal::is_in_mask<J> (exact(mask))));
  }

  template <typename I, typename J>
  sub_image<const I, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(const Image<I>&& ima, const Image<J>& mask)
  {
    static_assert(std::is_convertible<mln_value(J), bool>::value, "J's value type must be convertible to bool.");
    return make_subimage(move_exact(ima), rng::filter(exact(mask).domain(), internal::is_in_mask<J> (exact(mask))));
  }

  template <typename I, typename J>
  sub_image<I, filtered_range<const typename J::domain_type&, internal::is_in_mask<J> > >
  make_subimage(Image<I>&& ima, const Image<J>& mask)
  {
    static_assert(std::is_convertible<mln_value(J), bool>::value, "J's value type must be convertible to bool.");
    return make_subimage(move_exact(ima), rng::filter(exact(mask).domain(), internal::is_in_mask<J> (exact(mask))));
  }


  // Operator |
  template <typename I, typename DomainOrMask>
  auto operator| (const Image<I>& ima, const DomainOrMask& other)
    -> decltype( make_subimage(exact(ima), other) )
  {
    return make_subimage(exact(ima), other);
  }

  template <typename I, typename DomainOrMask>
  auto operator| (Image<I>& ima, const DomainOrMask& other)
    -> decltype( make_subimage(exact(ima), other) )
  {
    return make_subimage(exact(ima), other);
  }

  template <typename I, typename DomainOrMask>
  auto operator| (const Image<I>&& ima, const DomainOrMask& other)
    -> decltype( make_subimage(move_exact(ima), other) )
  {
    return make_subimage(move_exact(ima), other);
  }

  template <typename I, typename DomainOrMask>
  auto operator| (Image<I>&& ima, const DomainOrMask& other)
    -> decltype( make_subimage(move_exact(ima), other) )
  {
    return make_subimage(move_exact(ima), other);
  }



  /******************************************/
  /****              Traits              ****/
  /******************************************/


  template <typename Image, typename Domain>
  struct image_traits< sub_image<Image, Domain> >
  {
    typedef typename std::remove_reference<Image>::type I;

    typedef std::true_type      accessible;
    typedef forward_image_tag   category; // FIXME: category depends on domain category
    typedef std::integral_constant<bool,
                                   not std::is_reference<Image>::value and
                                   image_traits<I>::concrete::value> concrete;
    typedef std::false_type	indexable; // FIXME: depends
    typedef mln::extension::none_extension_tag extension;
  };


  template <typename I, typename Domain>
  struct image_concrete< sub_image<I, Domain> >
  {
    typedef sub_image<mln_concrete(I), Domain> type;
  };

  template <typename I, typename Domain, typename V>
  struct image_ch_value< sub_image<I, Domain>, V >
  {
    typedef sub_image<mln_ch_value(I, V), Domain> type;
  };

  // namespace internal
  // {
  //   template <class I, class Domain>
  //   struct image_init_from< sub_image<I, Domain> >
  //   {
  //     typedef sub_image<
  //       typename image_init_from<typename std::remove_reference<I>::type>::type,
  //       Domain> type;
  //   };

  // }

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename I, typename Domain>
  struct sub_image : image_base< sub_image<I, Domain>,
                                 typename range_value<Domain>::type,
                                 typename image_value<typename std::remove_reference<I>::type>::type >
  {
    BOOST_CONCEPT_ASSERT((AccessibleImage<typename std::decay<I>::type>));

  private:
    typedef typename std::remove_reference<I>::type  image_t;
    typedef sub_image<I, Domain>		     self_t;

    static_assert( std::is_convertible<typename range_value<Domain>::type, typename image_t::point_type>::value,
		"Domain's site type must be convertible to image's site type." );


    // static_assert( std::is_convertible<
    //                std::reference_wrapper<image_t>,
    //                std::reference_wrapper<const image_t> >::value, "pas conv");


    typedef transform_iterator<typename Domain::iterator,
                               std::reference_wrapper<image_t> >           value_iterator;
    typedef transform_iterator<typename Domain::iterator,
                               std::reference_wrapper<const image_t> >     const_value_iterator;

    struct pix_fun_t
    {
      pix_fun_t(image_t& ima)
        : m_ima(&ima)
      {
      }

      typename image_pixel<image_t>::type
      operator() (const mln_point(I)& p) const
      {
        return m_ima->pixel(p);
      }
    private:
      friend struct sub_image::const_pix_fun_t;
      image_t* m_ima;
    };

    struct const_pix_fun_t
    {
      const_pix_fun_t(const image_t& ima)
        : m_ima(&ima)
      {
      }

      const_pix_fun_t(const pix_fun_t& other)
        : m_ima(other.m_ima)
      {
      }

      typename image_const_pixel<image_t>::type
      operator() (const mln_point(I)& p) const
      {
        return m_ima->pixel(p);
      }
    private:
      const image_t* m_ima;
    };

    typedef rebind_pixel_iterator<
      self_t, transform_iterator<
                typename Domain::iterator,
                pix_fun_t > >                     pixel_iterator;
    typedef rebind_pixel_iterator<
      const self_t, transform_iterator<
                      typename Domain::iterator,
                      const_pix_fun_t > >         const_pixel_iterator;


  public:
    typedef typename range_value<Domain>::type                  point_type;
    typedef typename image_value<image_t>::type                 value_type;
    typedef typename image_reference<image_t>::type             reference;
    typedef typename image_const_reference<image_t>::type       const_reference;
    typedef Domain						domain_type;

    typedef rebinded_pixel<self_t, typename image_pixel<image_t>::type>               pixel_type;
    typedef rebinded_pixel<const self_t, typename image_const_pixel<image_t>::type>   const_pixel_type;

    typedef iterator_range<value_iterator>		value_range;
    typedef iterator_range<const_value_iterator>	const_value_range;
    typedef iterator_range<pixel_iterator>		pixel_range;
    typedef iterator_range<const_pixel_iterator>	const_pixel_range;

    template <class, class>
    friend struct sub_image;

    sub_image(I&& ima, const Domain& domain)
      : m_ima(std::forward<I>(ima)), m_domain(domain)
    {
    }

    sub_image() = default;


    template <typename OtherImage, typename OtherDomain>
    sub_image(const sub_image<OtherImage, OtherDomain>& other, mln::init)
      : m_ima(imchvalue<value_type>(other.m_ima)),
        m_domain(other.m_domain)
    {
    }

    template <typename OtherImage, typename OtherDomain>
    sub_image(const sub_image<OtherImage, OtherDomain>& other, const value_type& v)
      : m_ima(imchvalue<value_type>(other.m_ima).init(v)),
        m_domain(other.m_domain)
    {
    }

    const Domain& domain() const
    {
      return m_domain;
    }

    value_range values()
    {
      //return make_iterator_range( value_iterator(m_ima, m_domain.iter()) );
      return value_range(value_iterator(m_domain.iter(), m_ima));

    }

    const_value_range values() const
    {
      //return make_iterator_range( const_value_iterator(m_ima, m_domain.iter()) );
      return const_value_range(const_value_iterator(m_domain.iter(), m_ima));
    }

    pixel_range pixels()
    {
      //return make_iterator_range( pixel_iterator(m_ima, m_domain.iter(), *this) );
      return pixel_range(pixel_iterator
                         (*this, make_transform_iterator
                          (m_domain.iter(),
                           pix_fun_t(m_ima))));
    }

    const_pixel_range pixels() const
    {
      return const_pixel_range(const_pixel_iterator
                               (*this, make_transform_iterator
                                (m_domain.iter(), const_pix_fun_t(m_ima))));
    }

    reference
    operator() (const point_type& p)
    {
      mln_precondition(rng::has(m_domain, p));
      mln_precondition(rng::has(m_ima.domain(), p));
      return m_ima(p);
    }

    const_reference
    operator() (const point_type& p) const
    {
      mln_precondition(rng::has(m_domain, p));
      mln_precondition(rng::has(m_ima.domain(), p));
      return m_ima(p);
    }

    reference
    at (const point_type& p)
    {
      return m_ima.at(p);
    }

    const_reference
    at (const point_type& p) const
    {
      return m_ima.at(p);
    }

    pixel_type
    pixel_at (const point_type& p)
    {
      return pixel_type(*this, m_ima.pixel_at(p));
    }

    const_pixel_type
    pixel_at (const point_type& p) const
    {
      return const_pixel_type(*this, m_ima.pixel_at(p));
    }

    pixel_type
    pixel (const point_type& p)
    {
      mln_precondition(rng::has(m_domain, p));
      return pixel_type(*this, m_ima.pixel_at(p));
    }

    const_pixel_type
    pixel (const point_type& p) const
    {
      mln_precondition(rng::has(m_domain, p));
      return const_pixel_type(*this, m_ima.pixel_at(p));
    }


  private:
    I           m_ima;
    Domain      m_domain;
  };


} // end of namespace mln

# include <mln/core/image/sub_image.spe.hpp>

#endif //!MLN_CORE_IMAGE_SUB_IMAGE_HPP
