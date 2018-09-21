#ifndef MLN_CORE_IMAGE_IMAGE_HPP
# warning "Should not been included as a standalone"
# include <mln/core/image/image.hpp>
#endif

#ifndef MLN_CORE_IMAGE_INTERNAL_INITIALIZER_HPP
# define MLN_CORE_IMAGE_INTERNAL_INITIALIZER_HPP

# include <mln/core/image/internal/reindex.hpp>


namespace mln
{

  namespace internal
  {

    template <class OutputImage, class InputImageOrDomain>
    struct initializer;

    // Traits
    template <class I>
    struct image_init_from
    {
      typedef I type;
    };
  }

  struct init {};

  template <typename V, class I>
  internal::initializer<mln_ch_value(I, V), I>
  imchvalue(const Image<I>& ref);

  template <class I>
  internal::initializer<mln_concrete(I), I>
  imconcretize(const Image<I>& ref);


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {
    /// FWD declaration
    template <typename Neighborhood>
    int
    get_border_from_nbh(const Neighborhood& nbh);


    template <class OutputImage, class InputImageOrDomain>
    struct initializer
    {
    private:
      static constexpr bool has_border = image_has_border<OutputImage>::value;

    public:
      initializer(const InputImageOrDomain& ref)
        : m_ref(ref),
          m_status(nullptr),
          m_has_init(false)
      {
        __init_border();
      }

      initializer(initializer&&) = default;
      initializer& operator=(initializer&&) = default;

      initializer&
      init(const mln_value(OutputImage)& v)
      {
        m_has_init = true;
        m_init = v;
        return *this;
      }

      template <class dummy=initializer&>
      typename std::enable_if<has_border, dummy>::type
      border(unsigned b)
      {
        if (b > m_border)
          m_border = b;
        return *this;
      }

      template <class N, class dummy=initializer&>
      typename std::enable_if<has_border, dummy>::type
      adjust(const Neighborhood<N>& nbh)
      {
        int b = mln::internal::get_border_from_nbh(exact(nbh));
        if (b >= 0) {
          if (b > m_border)
            m_border = b;
        } else if (m_status)
          *m_status = -1;
        return *this;
      }

      operator OutputImage() const
      {
        OutputImage out = this->do_it();
        this->__reindex(out);
        return out;
      }

      initializer&
      get_status(int& status)
      {
        m_status = &status;
        *m_status = 0;
        return *this;
      }

    private:
      // prevent auto ima = imconcrete(...)
      initializer(const initializer&) = delete;
      initializer& operator=(const initializer&) = delete;

    protected:
      static constexpr bool is_input_an_image = mln::is_a<InputImageOrDomain, mln::Image>::value;
      static constexpr bool input_has_border =
        std::conditional< is_input_an_image,
                          image_has_border<InputImageOrDomain>,
                          std::false_type >::type::value;

      template <class dummy=void>
      OutputImage
      do_it(typename std::enable_if<( (dummy*)0 or (has_border and is_input_an_image)), dummy>::type* = nullptr) const
      {
        if (m_border != -1)
          {
            if (m_has_init)
              return OutputImage(m_ref, m_border, m_init);
            else
              return OutputImage(m_ref, m_border);
          }
        else
          {
            if (m_has_init)
              return OutputImage(m_ref, 3, m_init);
            else
              return OutputImage(m_ref, mln::init());
          }
      }

      template <class dummy=void>
      OutputImage
      do_it(typename std::enable_if<( (dummy*)0 or (has_border and !is_input_an_image)), dummy>::type* = nullptr) const
      {
        if (m_border != -1)
          {
            if (m_has_init)
              return make_image_from_domain<mln_value(OutputImage)> (m_ref, m_border, m_init);
            else
              return make_image_from_domain<mln_value(OutputImage)>(m_ref, m_border);
          }
        else
          {
            if (m_has_init)
              return make_image_from_domain<mln_value(OutputImage)>(m_ref, 3, m_init);
            else
              return make_image_from_domain<mln_value(OutputImage)>(m_ref);
          }
      }

      template <class dummy=void>
      OutputImage
      do_it(typename std::enable_if<( (dummy*)0 or (!has_border and is_input_an_image)), dummy>::type* = nullptr) const
      {
        if (m_has_init)
          return OutputImage(m_ref, m_init);
        else
          return OutputImage(m_ref, mln::init());
      }

      template <class dummy=void>
      OutputImage
      do_it(typename std::enable_if<( (dummy*)0 or (!has_border and !is_input_an_image)), dummy>::type* = nullptr) const
      {
        if (m_has_init)
          return make_image_from_domain<mln_value(OutputImage)>(m_ref, m_init);
        else
          return make_image_from_domain<mln_value(OutputImage)>(m_ref);
      }


      template <class dummy = void>
      typename std::enable_if<input_has_border,dummy>::type
      __init_border()
      {
        this->m_border = this->m_ref.border();
      }

      template <class dummy = void>
      typename std::enable_if<not input_has_border,dummy>::type
      __init_border()
      {
        this->m_border = -1;
      }

      template <class dummy = void>
      typename std::enable_if<is_input_an_image, dummy>::type
      __reindex(OutputImage& f) const
      {
        reindex(f, m_ref);
      }

      template <class dummy = void>
      typename std::enable_if<!is_input_an_image, dummy>::type
      __reindex(OutputImage&) const
      {
      }

    protected:
      InputImageOrDomain      m_ref;
      int*                    m_status;
      bool                    m_has_init;
      mln_value(OutputImage)  m_init;
      int                     m_border;
    };

  } // end of namespace mln::internal

  template <typename V, class I>
  internal::initializer<mln_ch_value(I, V), I>
  imchvalue(const Image<I>& _ref)
  {
    const I& ref = exact(_ref);
    return internal::initializer<mln_ch_value(I, V), I>(ref);
  }

  template <class I>
  internal::initializer<mln_concrete(I), I>
  imconcretize(const Image<I>& _ref)
  {
    const I& ref = exact(_ref);
    return internal::initializer<mln_concrete(I), I>(ref);
  }

} // end of namespace mln

#endif //!MLN_CORE_IMAGE_INTERNAL_INITIALIZER_HPP
