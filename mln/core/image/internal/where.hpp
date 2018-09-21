#ifndef MLN_CORE_IMAGE_INTERNAL_WHERE_HPP
# define MLN_CORE_IMAGE_INTERNAL_WHERE_HPP

# include <mln/core/iref.hpp>
# include <mln/core/pixel_utility.hpp>
# include <mln/core/iterator/filter_iterator.hpp>
# include <mln/core/iterator/transform_iterator.hpp>

namespace mln
{

  namespace internal
  {
    // FWD declaration
    template <typename I, class Predicate, class use_pix =
              typename internal::use_pix_helper<typename std::remove_reference<I>::type, Predicate>::type >
    struct where_t;

    template <typename I>
    struct where_binary_t;
  }

  /// \brief Return the domain of the image where the predicate is true.
  ///
  /// \param ima The input image E → V
  ///
  /// \param pred A unary predicate that takes either the pixel or the pixel
  /// value as argument i.e V → T or (E,V) → T, where \p T is convertible to
  /// bool.
  ///
  /// \return The subdomain of the image where the predicate is true. The result
  /// type fits the Domain concept.
  ///
  template <typename I, class Predicate>
  internal::where_t<const I&, Predicate>
  where(const Image<I>& ima, const Predicate& pred);

  template <typename I, class Predicate>
  internal::where_t<I, Predicate>
  where(Image<I>&& ima, const Predicate& pred);


  /// \brief  Return the domain of the image where it is true.
  ///
  /// \param ima The input image
  ///
  /// \return The subdomain of the image where it is true. The result
  /// type fits the Domain concept.
  ///
  /// \pre The value type of the image must be convertible to bool.
  template <typename I>
  internal::where_binary_t<const I&>
  where(const Image<I>& ima);

  template <typename I>
  internal::where_binary_t<I>
  where(Image<I>&& ima);


  /******************************************/
  /**** Specialization image from domain ****/
  /******************************************/

  // FWD declaration
  template <class I, class Domain>
  struct sub_image;

  template <class I>
  struct image_from_domain< internal::where_binary_t<I> >
  {
    template <typename V>
    struct apply
    {
      typedef sub_image< mln_ch_value(I,V), internal::where_binary_t<I> > type;
    };
  };

  template <class I, class Predicate>
  struct image_from_domain< internal::where_t<I, Predicate> >
  {
    template <typename V>
    struct apply
    {
      typedef sub_image< mln_ch_value(I,V), internal::where_t<I, Predicate> > type;
    };
  };

  template <class V, class I, class Predicate>
  sub_image< mln_ch_value(I,V), internal::where_t<I, Predicate> >
  make_image_from_domain(const internal::where_t<I, Predicate>& domain)
  {
    typedef sub_image< mln_ch_value(I,V), internal::where_t<I, Predicate> > OutputImage;
    return OutputImage(imchvalue<V>(domain.image()), domain);
  }

  template <class V, class I, class Predicate>
  sub_image< mln_ch_value(I,V), internal::where_t<I, Predicate> >
  make_image_from_domain(const internal::where_t<I, Predicate>& domain, const V& v)
  {
    typedef sub_image< mln_ch_value(I,V), internal::where_t<I, Predicate> > OutputImage;
    return OutputImage(imchvalue<V>(domain.image()).init(v), domain);
  }

  template <class V, class I>
  sub_image< mln_ch_value(I,V), internal::where_binary_t<I> >
  make_image_from_domain(const internal::where_binary_t<I>& domain)
  {
    typedef sub_image< mln_ch_value(I,V), internal::where_binary_t<I> > OutputImage;
    return OutputImage(imchvalue<V>(domain.image()), domain);
  }

  template <class V, class I>
  sub_image< mln_ch_value(I,V), internal::where_binary_t<I> >
  make_image_from_domain(const internal::where_binary_t<I>& domain, const V& v)
  {
    typedef sub_image< mln_ch_value(I,V), internal::where_binary_t<I> > OutputImage;
    return OutputImage(imchvalue<V>(domain.image()).init(v), domain);
  }


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {

    // Specialization when the Predicate uses pixels as argument.
    template <typename I, typename Predicate>
    struct where_t<I, Predicate, std::true_type>
    {
      template <class, class, class>
      friend struct where_t;

      typedef transform_iterator<filter_iterator<mln_cpxter(I), Predicate>,
                                 internal::meta_get_pixel_point>        iterator;
      typedef iterator                                                  const_iterator;

      where_t() = default;
      where_t(const where_t&) = default;
      where_t(where_t&&) = default;

      where_t(I&& ima, const Predicate& pred)
        : m_ima(std::forward<I>(ima)),
          m_pred(pred)
      {
      }

      template <class J, class Pred2>
      where_t(const where_t<J, Pred2>& other,
              typename std::enable_if<(std::is_convertible<J,I>::value and
                                       std::is_convertible<Pred2,Predicate>::value
                                       )>::type* = NULL)
        : m_ima (other.m_ima),
          m_pred (other.m_pred)
      {
      }


      iterator iter() const
      {
        return iterator(filter_iterator<mln_cpxter(I), Predicate>(m_ima.get().pixels().iter(), m_pred),
                        internal::meta_get_pixel_point());
      }

      bool has(const mln_point(I)& p) const
      {
        return m_ima.get().domain().has(p) and m_pred(m_ima.get().pixel(p));
      }

      const I& image() const
      {
        return m_ima.get();
      }

    private:
      mln::iref<I&&>    m_ima;
      Predicate         m_pred;
    };


    // Specialization when the Predicate uses values as argument.
    template <typename I, typename Predicate>
    struct where_t<I, Predicate, std::false_type>
    {
      template <class, class, class>
      friend struct where_t;

    private:
      struct predicate_t
      {
        bool
        operator() (const mln_cpixel(I)& px) const
        {
          return pred(px.val());
        }

        Predicate pred;
      };

    public:
      typedef transform_iterator<filter_iterator<mln_cpxter(I), predicate_t>,
                                 internal::meta_get_pixel_point>        iterator;
      typedef iterator                                                  const_iterator;

      where_t() = default;
      where_t(const where_t&) = default;
      where_t(where_t&&) = default;

      where_t(I&& ima, const Predicate& pred)
        : m_ima(std::forward<I>(ima)),
          m_pred{pred}
      {
      }

      template <class J, class Pred2>
      where_t(const where_t<J, Pred2>& other,
              typename std::enable_if<(std::is_convertible<J,I>::value and
                                       std::is_convertible<Pred2,Predicate>::value
                                       )>::type* = NULL)
        : m_ima (other.m_ima),
          m_pred (other.m_pred)
      {
      }

      iterator iter() const
      {
        return iterator(filter_iterator<mln_cpxter(I), predicate_t>(m_ima.get().pixels().iter(), m_pred),
                        internal::meta_get_pixel_point());
      }

      bool has(const mln_point(I)& p) const
      {
        return m_ima.get().domain().has(p) and m_pred.pred(m_ima.get()(p));
      }

      const I& image() const
      {
        return m_ima.get();
      }

    private:
      mln::iref<I&&>      m_ima;
      predicate_t         m_pred;
    };


    template <typename I>
    struct where_binary_t
    {
    public:
      typedef transform_iterator<filter_iterator<mln_cpxter(I), meta_get_pixel_value>,
                                 meta_get_pixel_point>          iterator;
      typedef iterator                                          const_iterator;

      where_binary_t() = default;

      where_binary_t(I&& ima)
        : m_ima(std::forward<I>(ima))
      {
      }

      iterator iter() const
      {
        return iterator(filter_iterator<mln_cpxter(I), meta_get_pixel_value>(m_ima.get().pixels().iter(),
                                                                             meta_get_pixel_value()),
                         internal::meta_get_pixel_point() );
      }

      bool has(const mln_point(I)& p) const
      {
        return m_ima.get().domain().has(p) and m_ima.get()(p);
      }

      const I& image() const
      {
        return m_ima.get();
      }

    private:
      mln::iref<I&&>         m_ima;
    };

  }

  template <typename I, class Predicate>
  internal::where_t<const I&, Predicate>
  where(const Image<I>& ima, const Predicate& pred)
  {
    return internal::where_t<const I&, Predicate>(exact(ima), pred);
  }

  template <typename I>
  internal::where_binary_t<const I&>
  where(const Image<I>& ima)
  {
    return internal::where_binary_t<const I&>(exact(ima));
  }

  template <typename I, class Predicate>
  internal::where_t<I, Predicate>
  where(Image<I>&& ima, const Predicate& pred)
  {
    return internal::where_t<I, Predicate>(move_exact(ima), pred);
  }

  template <typename I>
  internal::where_binary_t<I>
  where(Image<I>&& ima)
  {
    return internal::where_binary_t<I>(move_exact(ima));
  }

} // end of namespace mln

#endif //!MLN_CORE_IMAGE_INTERNAL_WHERE_HPP
