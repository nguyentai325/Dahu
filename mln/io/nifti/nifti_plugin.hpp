#ifndef MLN_IO_NIFTI_NIFTI_PLUGIN_HPP
# define MLN_IO_NIFTI_NIFTI_PLUGIN_HPP

# include <mln/io/plugin.hpp>
# include <mln/io/stream_wrapper.hpp>
# include <mln/core/grays.hpp>
# include <mln/core/colors.hpp>
# include <nifti1_io.h>

namespace mln
{

  namespace io
  {

    class nifti_reader_plugin : public PluginReader
    {
    public:
      nifti_reader_plugin();
      ~nifti_reader_plugin();

      virtual bool            support_istream() const final { return false; }

      virtual void            open(std::istream&) final;
      virtual void            open(const char* filename) final;
      virtual std::type_index get_value_type_id() const final;
      virtual int             get_bpp() const final;
      virtual box2d           get_domain() const final;

      virtual std::function<void(void*)> get_read_next_pixel_method() const final;
      virtual std::function<void(void*)> get_read_next_line_method() const final;

    private:
      void read_next_line_(void* out);
      void read_next_pixel_(void* out);


    protected:
      nifti_image*      m_nim;
      int               m_bpp;
      box2d             m_domain;
      std::type_index   m_vtype;

      std::function<void(void*)> m_read_next_line;
      std::function<void(void*)> m_read_next_pixel;
      char*              m_ptr;
      int                x,y;
      unsigned           m_stride;
    };


    /******************************/
    /***  Implementation        ***/
    /******************************/

    inline
    nifti_reader_plugin::nifti_reader_plugin()
      : m_vtype(typeid(void))
    {
    }

    nifti_reader_plugin::~nifti_reader_plugin()
    {
      if (m_nim)
        nifti_image_free(m_nim);
    }

    inline
    void
    nifti_reader_plugin::open(std::istream&)
    {
      throw MLNIOException("NIFTI does not support istreams as input.");
    }

    inline
    void
    nifti_reader_plugin::open(const char* filename)
    {
      // Only read header
      m_nim = nifti_image_read(filename, 0);

      if (not m_nim)
        throw MLNIOException("Unable to read the file");

      if (m_nim->ndim != 2) {
        std::stringstream msg;
        msg << "Only 2D images supported (for now) and it is " << m_nim->ndim << "D";
        msg << "Only the first two dimension will be loaded.";
        std::cerr << msg.str() << std::endl;
        //throw MLNIOException(msg.str());
      }

      m_bpp = m_nim->nbyper;

      switch (m_nim->datatype) {
        case DT_UINT8:  m_vtype = typeid(uint8);        break;
        case DT_UINT16: m_vtype = typeid(uint16);       break;
        case DT_UINT32: m_vtype = typeid(uint32);       break;
        case DT_INT8:  m_vtype = typeid(int8);        break;
        case DT_INT16: m_vtype = typeid(int16);       break;
        case DT_INT32: m_vtype = typeid(int32);       break;
        case DT_FLOAT32: m_vtype = typeid(float);       break;
        case DT_FLOAT64: m_vtype = typeid(double);       break;
        default: goto error;
      }

      m_read_next_line = std::bind(&nifti_reader_plugin::read_next_line_,
                                   this, std::placeholders::_1);
      m_read_next_pixel = std::bind(&nifti_reader_plugin::read_next_pixel_,
                                    this, std::placeholders::_1);


      m_domain.pmin = {0, 0};
      m_domain.pmax = {m_nim->ny, m_nim->nx};
      x = 0;

      nifti_image_load(m_nim);

      // FIXME: Should we go to the end
      // Maybe there is a bit for the orientation
      m_stride = m_bpp * m_domain.pmax[1];
      m_ptr = (char*)m_nim->data;
      m_ptr += (m_domain.pmax[0] - 1) * m_stride;

      return;

    error:
      std::stringstream msg;
      msg << "Unable to handle NIFTI file with value type: " << m_nim->datatype << "\n"
             << "bpp: " << m_bpp << "\n";
      throw MLNIOException(msg.str());
    }

    inline
    box2d
    nifti_reader_plugin::get_domain() const
    {
      return m_domain;
    }

    inline
    int
    nifti_reader_plugin::get_bpp() const
    {
      return m_bpp;
    }


    inline
    std::type_index
    nifti_reader_plugin::get_value_type_id() const
    {
      return m_vtype;
    }

    inline
    std::function<void(void*)>
    nifti_reader_plugin::get_read_next_line_method() const
    {
      return m_read_next_line;
    }

    inline
    std::function<void(void*)>
    nifti_reader_plugin::get_read_next_pixel_method() const
    {
      return m_read_next_pixel;
    }


    inline
    void nifti_reader_plugin::read_next_line_(void* out)
    {
      std::memcpy(out, m_ptr, m_domain.pmax[1] * m_bpp);
      m_ptr -= m_stride;
    }

    inline
    void nifti_reader_plugin::read_next_pixel_(void* out)
    {
      std::memcpy(out, m_ptr + x * m_bpp, m_bpp);
      ++x;
      if (x == m_domain.pmax[1]) {
        m_ptr -= m_stride;
        x = 0;
      }
    }

  }

}

#endif // ! MLN_IO_NIFTI_NIFTI_PLUGIN_HPP
