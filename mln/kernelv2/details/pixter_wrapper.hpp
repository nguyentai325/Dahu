#ifndef MLN_KERNELV2_DETAILS_PIXTER_WRAPPER_HPP
# define MLN_KERNELV2_DETAILS_PIXTER_WRAPPER_HPP

namespace mln
{
  namespace kernel
  {
    namespace details
    {

      template <typename I, typename Pixel>
      struct wrap_pixel :
        morpher_pixel_base< wrap_pixel<I, Pixel>,
                            typename std::remove_reference<Pixel>::type, I>
      {
        friend struct mln::morpher_core_access;
        typedef I                      image_type;
        typedef typename I::reference  reference;
        typedef typename I::value_type value_type;

        wrap_pixel(I& ima, const Pixel& pix)
          : m_ima(&ima), m_pix(pix)
        {
        }

        reference   val()   const { return m_pix.val(); }
        I& image() const          { return *m_ima; }

      private:
        I*    m_ima;
        Pixel m_pix;
      };



      template <class I, class Pixter>
      struct wrap_pixter :
        iterator_base< wrap_pixter<I, Pixter>,
                       wrap_pixel<I, typename Pixter::reference>,
                       wrap_pixel<I, typename Pixter::value_type> >
      {
        wrap_pixter(I& ima, Pixter* pixter)
          : m_ima(&ima), m_pixter(pixter)
        {
        }

        void init() { m_pixter->init(); }
        void next() { m_pixter->next(); }
        bool finished() const { return m_pixter->finished(); }

        wrap_pixel<I, typename Pixter::reference>
        dereference() const
        {
          return wrap_pixel<I, typename Pixter::reference>(*m_ima, *(*m_pixter));
        }

      private:
        I*      m_ima;
        Pixter* m_pixter;
      };

    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!MLN_KERNELV2_DETAILS_PIXTER_WRAPPER_HPP
