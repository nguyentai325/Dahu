#ifndef MLN_IO_LOADER_HPP
# define MLN_IO_LOADER_HPP

# include <type_traits>
# include <memory>
# include <algorithm>
# include <fstream>

# include <mln/core/value/value_traits.hpp>
# include <mln/core/image/image.hpp>
# include <boost/preprocessor/seq/for_each.hpp>
# include <mln/io/plugin.hpp>
# include <mln/io/ioexception.hpp>
# include <mln/io/internal/demangle.hpp>

namespace mln
{

  namespace io
  {

    template <class I>
    class Loader
    {
      BOOST_CONCEPT_ASSERT((Image<I>));

    public:
      void load(std::istream& s, Image<I>& ima, PluginReader* plugin,
                bool permissive);
      void load(const std::string& filename, Image<I>& ima, PluginReader* plugin,
                bool permissive);


    protected:
      void m_check_value_type_compatible() const;
      virtual bool is_value_type_convertible(std::type_index x,
                                             std::type_index y) const;
      virtual bool is_value_type_convertible_from(std::type_index src_type) const;

      virtual void m_resize() = 0;
      virtual void m_load();

      PluginReader* m_plugin;
      I*            m_ima;
      bool          m_permissive;
    };


    // Generic 2D loader
    template <class I, class category = typename image_traits<I>::category>
    class Loader2D : public Loader<I>
    {
      static_assert(std::is_same<typename I::domain_type, box2d>::value,
                    "The domain must be a box 2D.");

      virtual void m_resize();
    };

    // Specialization for fast images
    template <class I>
    class Loader2D<I, raw_image_tag> : public Loader<I>
    {
      static_assert(std::is_same<typename I::domain_type, box2d>::value,
                    "The domain must be a box 2D.");

      virtual void m_resize();
      virtual void m_load();
    };


    /***************************/
    /***  Implementation    ****/
    /***************************/

    template <class I>
    void
    Loader<I>::load(std::istream& s, Image<I>& ima, PluginReader* plugin,
                    bool permissive)
    {
      m_plugin = plugin;
      m_ima = & (exact(ima));
      m_permissive = permissive;

      m_plugin->initialize();
      m_plugin->open(s);
      this->m_check_value_type_compatible();
      this->m_resize();
      this->m_load();
    }

    template <class I>
    void
    Loader<I>::load(const std::string& s, Image<I>& ima, PluginReader* plugin,
                    bool permissive)
    {
      if (not plugin->support_istream() )
        {
          m_plugin = plugin;
          m_ima = & (exact(ima));
          m_permissive = permissive;

          m_plugin->initialize();
          m_plugin->open(s.c_str());

          this->m_check_value_type_compatible();
          this->m_resize();
          this->m_load();
        }
      else if (s == "-")
        {
          this->load(std::cin, ima, plugin, permissive);
        }
      else
        {
          std::ifstream stream(s, std::ios_base::binary);
          this->load(stream, ima, plugin, permissive);
        }
    }

    namespace internal
    {
      template <typename VIN, typename VOUT>
      void value_convert(void* buffer_in, void* buffer_out, std::size_t n)
      {
        VIN*  in = (VIN*) buffer_in;
        VOUT* out = (VOUT*) buffer_out;
        std::copy(in, in + n, out);
      }


# define MLN_INTERNAL_SCALAR_TYPE_SEQ (bool) (int8) (uint8) (int16) (uint16)\
        (int32) (uint32) (int64) (uint64) (float) (double)

# define MLN_INTERNAL_VEC_TYPE_SEQ (rgb8) (rgba8) (rgb16) (rgba16) (rgb32) (rgba32)


      template <typename VIN, typename VOUT, class Enable = void>
      struct default_converter
      {
        static constexpr std::nullptr_t ptr = nullptr;
      };

      template <typename VIN, typename VOUT>
      struct default_converter<VIN, VOUT, typename std::enable_if
                               <std::is_scalar<VIN>::value and std::is_scalar<VOUT>::value and
                                value_traits<VIN>::quant < value_traits<VOUT>::quant>::type>
      {
        static constexpr void (*ptr) (void*,void*,std::size_t) = &value_convert<VIN, VOUT>;
      };


# define MLN_INTERNAL_CONV_EXPAND(r, data, VIN)                         \
      if (sidx == typeid(VIN)) return default_converter<VIN, VOUT>::ptr;

      template <typename VOUT>
      std::function<void (void*, void*, std::size_t n)>
      defaut_value_type_converter(std::type_index sidx)
      {
        BOOST_PP_SEQ_FOR_EACH(MLN_INTERNAL_CONV_EXPAND, _, MLN_INTERNAL_SCALAR_TYPE_SEQ);

        return nullptr;
      }

    }

    template <class I>
    inline
    bool
    Loader<I>::is_value_type_convertible(std::type_index sidx,
                                         std::type_index tidx) const
    {

      return sidx == tidx;
    }

    template <class I>
    inline
    bool
    Loader<I>::is_value_type_convertible_from(std::type_index sidx) const
    {
      typedef mln_value(I) V;
      return sidx == typeid(V) or internal::defaut_value_type_converter<V>(sidx);
    }


    template <class I>
    void
    Loader<I>::m_check_value_type_compatible() const
    {
      typedef mln_value(I) V;
      std::type_index sidx = m_plugin->get_value_type_id();
      std::type_index tidx = typeid(V);

      if (sidx != tidx and (not m_permissive or
                            not is_value_type_convertible_from(sidx)))
        {
          std::string ex = "Value types incompatibles: ";
          (ex += "trying to load ") += internal::demangle(sidx.name());
          (ex += " in an image of ") += internal::demangle(tidx.name());
          throw MLNIOException(ex);
        }
    }

    template <class I>
    void
    Loader<I>::m_load()
    {
      std::type_index sidx = m_plugin->get_value_type_id();
      std::type_index tidx = typeid(mln_value(I));

      std::function<void (void*,void*,std::size_t)> from_to =
        internal::defaut_value_type_converter< mln_value(I) >(sidx);

      std::function<void (void*)> read_next_pixel =
        m_plugin->get_read_next_pixel_method();


      if (sidx != tidx and from_to) // need to convert
        {
          int bpp       = m_plugin->get_bpp();
          int bytes     = bpp / 8 + (bpp % 8 != 0);
          void* valptr  = std::malloc(bytes);

          mln_pixter(px, *m_ima);
          mln_forall(px)
          {
            read_next_pixel(valptr);
            from_to(valptr, &(px->val()), 1);
          }

          std::free(valptr);
        }
      else
        {
          mln_assertion(sidx == tidx);

          mln_value(I) v;
          mln_pixter(px, *m_ima);
          mln_forall(px)
          {
            read_next_pixel((void*)&v);
            px->val() = v;
          }
        }
      // else
      //   {
      //     std::string ex = "Value types incompatibles: ";
      //     (ex += "trying to load ") += internal::demangle(sidx.name());
      //     (ex += " in an image of ") += internal::demangle(tidx.name());
      //     throw MLNIOException(ex);
      //   }
    }


    template <class I, class category>
    void
    Loader2D<I, category>::m_resize()
    {
      PluginReader2D* plug = reinterpret_cast<PluginReader2D*>(this->m_plugin);
      box2d dom = plug->get_domain();
      this->m_ima->resize(dom);
    }

    template <class I>
    void
    Loader2D<I, raw_image_tag>::m_resize()
    {
      PluginReader2D* plug = reinterpret_cast<PluginReader2D*>(this->m_plugin);
      box2d dom = plug->get_domain();
      this->m_ima->resize(dom, this->m_ima->border());
    }

    template <class I>
    void
    Loader2D<I, raw_image_tag>::m_load()
    {
      PluginReader2D* plug = reinterpret_cast<PluginReader2D*>(this->m_plugin);
      box2d dom = this->m_ima->domain();
      point2d p = dom.pmin;
      point2d q = dom.pmax;

      std::function<void (void*)> read_next_line =
        plug->get_read_next_line_method();

      std::type_index sidx = plug->get_value_type_id();
      std::type_index tidx = typeid(mln_value(I));

      std::function<void (void*,void*,std::size_t)> from_to =
        internal::defaut_value_type_converter<mln_value(I)>(sidx);

      if (sidx != tidx and from_to) // we need to convert the value type
        {
          int bpp       = plug->get_bpp();
          std::size_t n = q[1] - p[1];
          void* buffer  = std::malloc(bpp * n);

          for (; p[0] != q[0]; ++p[0])
             {
               void* lineptr = &(*this->m_ima)(p);
               read_next_line(buffer);
               from_to(buffer, lineptr, n);
             }

          std::free(buffer);
        }
      else
        {
          mln_assertion(sidx == tidx);
          for (; p[0] != q[0]; ++p[0])
            {
              void* lineptr = &(*this->m_ima)(p);
              read_next_line(lineptr);
            }
        }
    }

  }

}

#endif // ! MLN_IO_LOADER_HPP
