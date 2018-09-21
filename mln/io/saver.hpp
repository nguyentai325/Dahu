#ifndef MLN_IO_SAVER_HPP
# define MLN_IO_SAVER_HPP

# include <mln/core/image/image.hpp>
# include <mln/io/plugin.hpp>
# include <mln/io/ioexception.hpp>
# include <mln/io/internal/demangle.hpp>


namespace mln
{
  namespace io
  {

    template <class I>
    class Saver
    {
      BOOST_CONCEPT_ASSERT((Image<I>));

    public:
      void save(const Image<I>& ima, PluginWriter* plugin, bool permissive);

    protected:
      virtual void m_set_domain(const I& ima, PluginWriter* plugin) = 0;
      virtual void m_save(const I& ima, PluginWriter* plugin, bool permissive);
    };

    template <class I, class category = typename image_traits<I>::category>
    class Saver2D : public Saver<I>
    {
      static_assert(std::is_same<typename I::domain_type, box2d>::value,
                    "The domain must be a box2d.");

    protected:
      virtual void m_set_domain(const I& ima, PluginWriter* plugin);
    };

    template <class I>
    class Saver2D<I, raw_image_tag> : public Saver<I>
    {
      static_assert(std::is_same<typename I::domain_type, box2d>::value,
                    "The domain must be a box2d.");

    protected:
      virtual void m_set_domain(const I& ima, PluginWriter* plugin);
      virtual void m_save(const I& ima, PluginWriter* plugin, bool permissive);
    };

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    template <class I>
    void
    Saver<I>::save(const Image<I>& ima_, PluginWriter* plugin,
                   bool permissive)
    {
      const I& ima = exact(ima_);

      if (not plugin->can_write(typeid(mln_value(I))))
        {
          std::string msg = "The plugin does not support writing " +
            internal::demangle(typeid(mln_value(I)).name());
          throw MLNIOException(msg);
        }

      plugin->set_value_type(typeid(mln_value(I)));
      this->m_set_domain(ima, plugin);
      this->m_save(ima, plugin, permissive);
    }

    template <class I>
    void
    Saver<I>::m_save(const I& ima, PluginWriter* plugin, bool permissive)
    {
      (void) permissive;
      std::function<void (void*)> write_next_pixel =
        plugin->get_write_next_pixel_method();

      mln_foreach(mln_value(I) v, ima.values())
        {
          write_next_pixel((void*) &v);
        }
    }

    template <class I, class category>
    void
    Saver2D<I, category>::m_set_domain(const I& ima, PluginWriter* plugin)
    {
      PluginWriter2D* plug = reinterpret_cast<PluginWriter2D*>(plugin);
      plug->set_domain(ima.domain());
    }

    template <class I>
    void
    Saver2D<I, raw_image_tag>::m_set_domain(const I& ima, PluginWriter* plugin)
    {
      PluginWriter2D* plug = reinterpret_cast<PluginWriter2D*>(plugin);
      plug->set_domain(ima.domain());
    }

    template <class I>
    void
    Saver2D<I, raw_image_tag>::m_save(const I& ima, PluginWriter* plugin,
                                      bool permissive)
    {
      (void) permissive;

      PluginWriter2D* plug = reinterpret_cast<PluginWriter2D*>(plugin);
      std::function<void(void*)> write_line =
        plug->get_write_next_line_method();

      point2d p = ima.domain().pmin;
      point2d pmax = ima.domain().pmax;

      for (; p[0] != pmax[0]; ++p[0])
        {
          void* lineptr = (void*) &(ima(p));
          write_line(lineptr);
        }
    }


  } // end of namespace mln::io
} // end of namespace mln

#endif //!MLN_IO_SAVER_HPP
