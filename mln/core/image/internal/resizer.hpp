#ifndef MLN_CORE_IMAGE_IMAGE_HPP
# warning "Should not been included as a standalone"
# include <mln/core/image/image.hpp>
#endif

#ifndef MLN_CORE_IMAGE_INTERNAL_RESIZER_HPP
# define MLN_CORE_IMAGE_INTERNAL_RESIZER_HPP

# include <mln/core/image/internal/reindex.hpp>

namespace mln
{

  namespace internal
  {
    template <class OutputImage,
              class InputImageOrInputDomain>
    struct resizer;
  }

  /// \brief Initialize an image from an other such that
  /// domains correspond. If \p ima and \p other are indexable
  /// indexes match as well.
  ///
  /// This is equivalent to `ima = I(other, mln::init())` but requires
  /// the image to be default constructible. If not `I` is not default-constructible
  /// you better use `I ima = imchvalue<mln_value(I)>(other)`.
  ///
  /// \param ima Image to be initialized.
  /// \param other Image to initialize with.
  ///
  template <typename I, typename J>
  internal::resizer<I,J>
  resize(Image<I>& ima, const Image<J>& other); //__attribute__((deprecated));


  namespace internal
  {
    /// FWD declaration
    template <typename Neighborhood>
    int
    get_border_from_nbh(const Neighborhood& nbh);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace impl
    {
      template <class OutputImage, class R>
      struct resizer_base;

      template <class OutputImage,
                class InputDomain>
      struct resizer_by_domain;

      template <class OutputImage,
                class InputImage>
      struct resizer_by_image;

      template <class OutputImage, class R>
      struct resizer_base
      {
      private:
        static constexpr bool has_border = image_has_border<OutputImage>::value;

      public:
        resizer_base(OutputImage& out, const R& ref)
          : m_output(out),
            m_ref(ref),
            m_has_hook(true),
            m_failed(false),
            m_has_init(false)
        {
        }

        resizer_base(resizer_base&& other)
          : m_output(other.m_output),
            m_ref(other.m_ref),
            m_has_hook(true),
            m_failed(other.m_failed),
            m_has_init(other.m_has_init),
            m_init(other.m_init),
            m_border(other.m_border)
        {
          other.m_has_hook = false;
        }

        ~resizer_base()
        {
          do_it();
        }

        resizer_base&
        init(const mln_value(OutputImage)& v)
        {
          m_has_init = true;
          m_init = v;
          return *this;
        }

        template <class dummy=resizer_base&>
        typename std::enable_if<has_border, dummy>::type
        border(int b)
        {
          if (b > m_border)
            m_border = b;
          return *this;
        }

        template <class N, class dummy=resizer_base&>
        typename std::enable_if<has_border, dummy>::type
        adjust(const Neighborhood<N>& nbh)
        {
          int b = mln::internal::get_border_from_nbh(exact(nbh));
          if (b >= 0) {
            if (b > m_border)
              m_border = b;
          } else
            m_failed = true;
          return *this;
        }

        operator bool()
        {
          do_it();
          return !m_failed;
        }


      protected:
        resizer_base(const resizer_base&) = delete;
        resizer_base& operator=(const resizer_base&) = delete;


        template <class dummy=void>
        typename std::enable_if<has_border, dummy>::type
        do_it()
        {
          if (not m_has_hook)
            return;

          if (m_border != -1)
            {
              if (m_has_init)
                m_output.resize(m_ref.domain(), m_border, m_init);
              else
                m_output.resize(m_ref.domain(), m_border);
            }
          else
            {
              if (m_has_init)
                m_output.resize(m_ref.domain(), 3, m_init);
              else
                m_output.resize(m_ref.domain());
            }
          reindex(m_output, m_ref);
          m_has_hook = false;
        }

        template <class dummy=void>
        typename std::enable_if<!has_border, dummy>::type
        do_it()
        {
          if (not m_has_hook)
            return;

          if (m_has_init)
            m_output.resize(m_ref.domain(), m_init);
          else
            m_output.resize(m_ref.domain());

          reindex(m_output, m_ref);
          m_has_hook = false;
        }

      protected:
        OutputImage&            m_output;
        const R&                m_ref;

        bool                    m_has_hook;
        bool                    m_failed;
        bool                    m_has_init;
        mln_value(OutputImage)  m_init;
        int                     m_border;
      };

      template <class OutputImage,
                class InputDomain>
      struct resizer_by_domain : resizer_base<OutputImage, InputDomain>
      {
        resizer_by_domain(OutputImage& out, const InputDomain& domain)
          : resizer_base<OutputImage, InputDomain>(out, domain)
        {
          this->m_border = -1;
        }
      };

      template <class OutputImage,
                class InputImage>
      struct resizer_by_image : resizer_base<OutputImage, InputImage>
      {
        resizer_by_image(OutputImage& out, const InputImage& ref)
          : resizer_base<OutputImage, InputImage>(out, ref)
        {
          __init_border();
        }

      protected:
        static constexpr bool has_border = image_has_border<InputImage>::value;

        template <class dummy = void>
        typename std::enable_if<has_border,dummy>::type
        __init_border()
        {
          this->m_border = this->m_ref.border();
        }

        template <class dummy = void>
        typename std::enable_if<not has_border,dummy>::type
        __init_border()
        {
          this->m_border = -1;
        }

      };

    }

    template <class OutputImage,
              class InputImageOrInputDomain>
    struct resizer :
      std::conditional< is_a<InputImageOrInputDomain, Image>::value,
                        impl::resizer_by_image<OutputImage, InputImageOrInputDomain>,
                        impl::resizer_by_domain<OutputImage, InputImageOrInputDomain> >::type
    {
    private:
      typedef typename
      std::conditional< is_a<InputImageOrInputDomain, Image>::value,
                        impl::resizer_by_image<OutputImage, InputImageOrInputDomain>,
                        impl::resizer_by_domain<OutputImage, InputImageOrInputDomain> >::type
      base;

    public:
      resizer(OutputImage& out, const InputImageOrInputDomain& ref)
        : base(out, ref)
      {
      }
    };

  } // end of namespace mln::internal

  template <typename I, typename J>
  inline
  internal::resizer<I, J>
  resize(Image<I>& ima, const Image<J>& other)
  {
    return {exact(ima), exact(other)};
  }


} // end of namespace mln

#endif //!MLN_CORE_IMAGE_INTERNAL_RESIZER_HPP
