#ifndef MLN_IO_FREEIMAGE_PLUGIN_HPP
# define MLN_IO_FREEIMAGE_PLUGIN_HPP

# include <mln/io/plugin.hpp>
# include <mln/io/stream_wrapper.hpp>
# include <mln/core/grays.hpp>
# include <mln/core/colors.hpp>
# include <mln/colors/rgba.hpp>
# include <unordered_map>
# include <FreeImage.h>

namespace mln
{

  namespace io
  {

    class freeimage_reader_plugin : public PluginReader2D
    {
    public:
      freeimage_reader_plugin();
      ~freeimage_reader_plugin();

      virtual void open(std::istream& s) final;
      //virtual bool support_streaming() const;
      virtual void read_next_pixel(void* out) final;
      virtual void read_next_line(void* out) final;
      virtual std::function<void(void*)> get_read_next_pixel_method() const final;
      virtual std::function<void(void*)> get_read_next_line_method() const final;

      virtual std::type_index get_value_type_id() const final;

      virtual box2d get_domain() const final;
      virtual int   get_bpp() const final;

    private:
      void _load();
      template <class colortype, class FI_type> void read_next_line_rgb(void* out);
      template <class colortype, class FI_type> void read_next_pixel_rgb(void* out);
      template <class colortype, class FI_type> void read_next_line_rgba(void* out);
      template <class colortype, class FI_type> void read_next_pixel_rgba(void* out);
      void read_next_line_rgb8(void* out);
      void read_next_pixel_rgb8(void* out);
      void read_next_line_rgba8(void* out);
      void read_next_pixel_rgba8(void* out);
      void read_next_line_gray(void* out);
      void read_next_pixel_gray(void* out);
      void read_next_line_gray_palette(void* out);
      void read_next_pixel_gray_palette(void* out);
      void read_next_line_bool(void* out);
      void read_next_pixel_bool(void* out);
      void read_next_line_bool_inv(void* out);
      void read_next_pixel_bool_inv(void* out);

    private:
      std::function<void (void*)>   m_read_next_pixel;
      std::function<void (void*)> m_read_next_line;


      box2d             m_domain;
      FREE_IMAGE_FORMAT m_fif;
      FIBITMAP*         m_dib;
      char*             m_ptr;
      std::type_index   m_vtype;
      RGBQUAD*          m_palette;
      unsigned  pitch;          // Size in byte of the scanline (padding inc.)
      unsigned  byteperline;    // Size in byte of the line (without padding)
      int       bpp;            // Size in byte of the pixel (1/4/8bits bitmap => 1 bpp, 16bits -> 2bpp...)
      int       x, y; // internally used by read/write to store current position
    };

    class freeimage_writer_plugin : public PluginWriter2D
    {
    public:
      freeimage_writer_plugin(std::ostream& os, FREE_IMAGE_FORMAT fif = FIF_TIFF);
      ~freeimage_writer_plugin();

      //virtual bool support_streaming() const;
      virtual void write_next_pixel(void* src) final;
      virtual void write_next_line(void* src) final;
      virtual std::function<void(void*)> get_write_next_pixel_method() const final;
      virtual std::function<void(void*)> get_write_next_line_method() const final;

      virtual bool can_write(std::type_index) const final;
      virtual void set_domain(const box2d& domain) final;
      virtual void set_value_type(std::type_index type) final;

    private:
      void _allocate();
      void _get_fit_and_bpp_from_type(FREE_IMAGE_TYPE& fit, unsigned& bpp, std::type_index vtype) const;
      void write_next_line_rgb8(void* out);
      void write_next_pixel_rgb8(void* out);
      void write_next_line_rgba8(void* out);
      void write_next_pixel_rgba8(void* out);
      void write_next_line_gray(void* out);
      void write_next_pixel_gray(void* out);
      void write_next_line_bool(void* out);
      void write_next_pixel_bool(void* out);

    private:
      std::function<void (void*)>   m_write_next_pixel;
      std::function<void (void*)>   m_write_next_line;

      internal::ostream_wrapper m_osw;
      unsigned                  m_nrows;
      unsigned                  m_ncols;
      FREE_IMAGE_FORMAT         m_fif;
      FIBITMAP*                 m_dib;
      std::type_index           m_vtype;
      char*                     m_ptr;
      unsigned          pitch;          // Size in byte of the scanline (padding inc.)
      unsigned          byteperline;    // Size in byte of the line (without padding)
      unsigned          psz;            // Size in byte of a pixel
      unsigned          bpp;            // Size in bit of the pixel
      unsigned          x, y; // internally used by read/write to store current position
    };


    /************************/
    /*** Implementation   ***/
    /************************/

    inline
    freeimage_reader_plugin::freeimage_reader_plugin()
      : m_dib (NULL), m_vtype(typeid(void)), x (0), y (0)
    {
    }

    inline
    box2d
    freeimage_reader_plugin::get_domain() const
    {
      return m_domain;
    }

    inline
    int
    freeimage_reader_plugin::get_bpp() const
    {
      return bpp;
    }

    inline
    std::function<void(void*)>
    freeimage_reader_plugin::get_read_next_line_method() const
    {
      return m_read_next_line;
    }

    inline
    std::function<void(void*)>
    freeimage_reader_plugin::get_read_next_pixel_method() const
    {
      return m_read_next_pixel;
    }


    inline
    void
    freeimage_reader_plugin::open(std::istream& is)
    {
      mln_precondition(m_dib == NULL);


      internal::istream_wrapper isw(is);
      fi_handle handle = (fi_handle) &isw;

      FreeImageIO       fio = {
        internal::istream_wrapper::read,
        NULL,
        internal::istream_wrapper::seek,
        internal::istream_wrapper::tell
      };

      m_fif = FreeImage_GetFileTypeFromHandle(&fio, handle);

      if (m_fif == FIF_UNKNOWN || !FreeImage_FIFSupportsReading(m_fif))
        throw MLNIOException("File format not supported");

      m_dib = FreeImage_LoadFromHandle(m_fif, &fio, handle, 0);
      if (not m_dib)
        throw MLNIOException("Unable to read the file");

      this->_load();
      mln_postcondition(m_dib != NULL);
    }

    inline
    void
    freeimage_reader_plugin::_load()
    {
      mln_precondition(m_dib != NULL);
      m_domain.pmin[0] = 0;
      m_domain.pmin[1] = 0;
      m_domain.pmax[0] = FreeImage_GetHeight(m_dib);
      m_domain.pmax[1] = FreeImage_GetWidth(m_dib);

      pitch = FreeImage_GetPitch(m_dib);
      m_ptr = (char*)FreeImage_GetBits(m_dib) + pitch * (m_domain.pmax[0]-1);

      FREE_IMAGE_TYPE type = FreeImage_GetImageType(m_dib);
      int bppp = FreeImage_GetBPP(m_dib);
      int colortype = FreeImage_GetColorType(m_dib);
      bpp = bppp / 8;
      m_palette = FreeImage_GetPalette(m_dib);

      switch (type)
        {
          case FIT_BITMAP:
            break;
          case FIT_UINT16:
            m_vtype = typeid(uint16);
            break;
          case FIT_INT16:
            m_vtype = typeid(int16);
            break;
          case FIT_UINT32:
            m_vtype = typeid(uint32);
            break;
          case FIT_INT32:
            m_vtype = typeid(int32);
            break;
          case FIT_FLOAT:
            m_vtype = typeid(float);
            break;
          case FIT_DOUBLE:
            m_vtype = typeid(double);
            break;
          case FIT_RGB16:
            m_vtype = typeid(rgb16);
            break;
          case FIT_RGBA16:
            m_vtype = typeid(colors::rgba16);
            break;
          case FIT_RGBF:
            m_vtype = typeid(rgb<float>);
            break;
          case FIT_RGBAF:
            m_vtype = typeid(colors::rgba<float>);
            break;
          default:
            goto error;
        }
      // For any non-bitmap or 8-bit bitmap type, we simply memcpy the pixel
      // because the underlying FreeImage value type is the same as our.
      if (type != FIT_BITMAP or bppp == 8)
        {
          m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_gray,
                                       this, std::placeholders::_1);
          m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_gray,
                                        this, std::placeholders::_1);
        }

      // We only need to take care of those special cases:
      // * 8-bit RGB[A] (because Freeimage value type ordering is OS dependant)
      // * 1-bit monochrome bitmap
      // * 4-bit bitmap
      if (type == FIT_BITMAP)
        {
          switch (bppp)
            {
            case 1:
              m_vtype = typeid(bool);
              if (colortype == FIC_MINISBLACK) {
                m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_bool,
                                             this, std::placeholders::_1);
                m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_bool,
                                              this, std::placeholders::_1);
              } else {
                m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_bool_inv,
                                             this, std::placeholders::_1);
                m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_bool_inv,
                                              this, std::placeholders::_1);
              }
              break;
            case 8:
              if (colortype == FIC_MINISBLACK) {
                m_vtype = typeid(uint8);
                break;
              } else if (colortype == FIC_PALETTE) {
                m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_gray_palette,
                                             this, std::placeholders::_1);
                m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_gray_palette,
                                              this, std::placeholders::_1);
                m_vtype = typeid(rgb8);
                break;
              }
              goto error;
            case 24:
              if (colortype == FIC_RGB) {
                m_vtype = typeid(rgb8);
                m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_rgb8,
                                             this, std::placeholders::_1);
                m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_rgb8,
                                              this, std::placeholders::_1);
                break;
              }
              goto error;
            case 32:
              if (colortype == FIC_RGBALPHA) {
                m_vtype = typeid(colors::rgba8);
                m_read_next_line = std::bind(&freeimage_reader_plugin::read_next_line_rgba8,
                                             this, std::placeholders::_1);
                m_read_next_pixel = std::bind(&freeimage_reader_plugin::read_next_pixel_rgba8,
                                              this, std::placeholders::_1);
                break;
              }

            default:
              goto error;
            }
        }
      return;

    error:
      m_vtype = typeid(void);
      std::string cstr_ctype;
      switch (colortype) {
      case FIC_MINISBLACK:
        cstr_ctype = "FIC_MINISBLACK";
        break;
      case FIC_MINISWHITE:
        cstr_ctype = "FIC_MINISWHITE";
        break;
      case FIC_PALETTE:
        cstr_ctype = "FIC_PALETTE";
        break;
      case FIC_RGB:
        cstr_ctype = "FIC_RGB";
        break;
      case FIC_RGBALPHA:
        cstr_ctype = "FIC_RGBALPHA";
        break;
      case FIC_CMYK:
        cstr_ctype = "FIC_CMYK";
        break;
      };

      std::stringstream msg;
      msg << "Unable to handle FORMAT: " << FreeImage_GetFormatFromFIF(m_fif) << "\n"
          << "with color type: " << cstr_ctype << "\n"
          << "bpp: " << std::to_string(bppp) << "\n";
      throw MLNIOException(msg.str());
    }


    inline
    std::type_index
    freeimage_reader_plugin::get_value_type_id() const
    {
      mln_precondition(m_dib != NULL);
      return m_vtype;
    }

    /******************************************/
    /****         Reading methods         ****/
    /******************************************/

    // 1-bit DIBs are stored using each bit as an index into the color table. The most
    // significant bit is the leftmost pixel.
    inline
    void freeimage_reader_plugin::read_next_line_bool(void* out_)
    {
      mln_precondition(m_dib != NULL);

      char* out = (char*)out_;
      for (int y = 0, z = 0; y < m_domain.pmax[1]; y += 8, z += 1)
        for (int b = 0; b < 8; ++b)
          out[y+b] = ((m_ptr[z] & (0x80 >> b)) != 0);

      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_bool(void* out_)
    {
      mln_precondition(m_dib != NULL);

      char* out = (char*)out_;
      *out = ((m_ptr[x] & (0x80 >> y)) != 0);

      if (++y == 8) {
        y = 0;
        if (++x == (int)byteperline)
          {
            x = 0;
            m_ptr -= pitch;
          }
      }
    }

    inline
    void freeimage_reader_plugin::read_next_line_bool_inv(void* out_)
    {
      mln_precondition(m_dib != NULL);

      char* out = (char*)out_;
      for (int y = 0, z = 0; y < m_domain.pmax[1]; y += 8, z += 1)
        for (int b = 0; b < 8; ++b)
          out[y+b] = ((m_ptr[z] & (0x80 >> b)) == 0);

      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_bool_inv(void* out_)
    {
      mln_precondition(m_dib != NULL);

      char* out = (char*)out_;
      *out = ((m_ptr[x] & (0x80 >> y)) == 0);

      if (++y == 8) {
        y = 0;
        if (++x == (int)byteperline)
          {
            x = 0;
            m_ptr -= pitch;
          }
      }
    }


    // 8-bit DIBs are the easiest to store because each byte is an
    // index into the color table.
    // Non standard image types such as short, long, float or double do not have a color
    // table. Pixels are stored in a similar way as 8-bit DIB
    inline
    void freeimage_reader_plugin::read_next_line_gray(void* out)
    {
      std::memcpy(out, m_ptr, m_domain.pmax[1] * bpp);
      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_gray(void* out)
    {
      std::memcpy(out, m_ptr + y * bpp, bpp);
      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }


    // BITMAP with a color lookup table
    // It may be 1,4 or 8 bits images
    inline
    void freeimage_reader_plugin::read_next_line_gray_palette(void* out)
    {
      rgb8* buffer = (rgb8*) out;
      uint8* ptr = (uint8*) m_ptr;
      for (int y = 0; y < m_domain.pmax[1]; ++y, ++ptr)
        {
          buffer[y][0] = m_palette[*ptr].rgbRed;
          buffer[y][1] = m_palette[*ptr].rgbGreen;
          buffer[y][2] = m_palette[*ptr].rgbBlue;
        }
      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_gray_palette(void* out)
    {
      rgb8* pxout = (rgb8*) out;
      uint8* ptr = (uint8*) m_ptr;
      (*pxout)[0] = m_palette[*ptr].rgbRed;
      (*pxout)[1] = m_palette[*ptr].rgbGreen;
      (*pxout)[2] = m_palette[*ptr].rgbBlue;

      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }


    // 24-bit DIBs have every 3 bytes representing a color, using the
    // same ordering as the RGBTRIPLE structure.
    template <class colortype, class FI_type>
    inline
    void freeimage_reader_plugin::read_next_line_rgb(void* out)
    {
      colortype* buffer = (colortype*) out;
      FI_type* ptr = (FI_type*) m_ptr;
      for (int y = 0; y < m_domain.pmax[1]; ++y, ++ptr)
        {
          buffer[y][0] = ptr->red;
          buffer[y][1] = ptr->green;
          buffer[y][2] = ptr->blue;
        }
      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_line_rgb8(void* out)
    {
      rgb8* buffer = (rgb8*) out;
      for (int y = 0; y < m_domain.pmax[1]; y++)
        {
          buffer[y][0] = m_ptr[y * bpp + FI_RGBA_RED];
          buffer[y][1] = m_ptr[y * bpp + FI_RGBA_GREEN];
          buffer[y][2] = m_ptr[y * bpp + FI_RGBA_BLUE];
        }
      m_ptr -= pitch;
    }

    template <class colortype, class FI_type>
    inline
    void freeimage_reader_plugin::read_next_pixel_rgb(void* out)
    {
      colortype* pixel = (colortype*) out;
      FI_type* ptr = (FI_type*)(m_ptr);
      (*pixel)[0] = ptr[y].red;
      (*pixel)[1] = ptr[y].green;
      (*pixel)[2] = ptr[y].blue;

      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_rgb8(void* out)
    {
      rgb8* pixel = (rgb8*) out;
      (*pixel)[0] = m_ptr[y * bpp + FI_RGBA_RED];
      (*pixel)[1] = m_ptr[y * bpp + FI_RGBA_GREEN];
      (*pixel)[2] = m_ptr[y * bpp + FI_RGBA_BLUE];
      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }

    // 32-bit DIBs have every 4 bytes representing a color, using the
    // same ordering as the RGBQUAD structure.
    template <class colortype, class FI_type>
    inline
    void freeimage_reader_plugin::read_next_line_rgba(void* out)
    {
      colortype* buffer = (colortype*) out;
      FI_type* ptr = (FI_type*) m_ptr;
      for (int y = 0; y < m_domain.pmax[1]; ++y, ++ptr)
        {
          buffer[y][0] = ptr->red;
          buffer[y][1] = ptr->green;
          buffer[y][2] = ptr->blue;
          buffer[y][3] = ptr->alpha;
        }
      m_ptr -= pitch;
    }

    inline
    void freeimage_reader_plugin::read_next_line_rgba8(void* out)
    {
      colors::rgba8* buffer = (colors::rgba8*) out;
      for (int y = 0; y < m_domain.pmax[1]; y++)
        {
          buffer[y][0] = m_ptr[y * bpp + FI_RGBA_RED];
          buffer[y][1] = m_ptr[y * bpp + FI_RGBA_GREEN];
          buffer[y][2] = m_ptr[y * bpp + FI_RGBA_BLUE];
          buffer[y][3] = m_ptr[y * bpp + FI_RGBA_ALPHA];
        }
      m_ptr -= pitch;
    }


    template <class colortype, class FI_type>
    inline
    void freeimage_reader_plugin::read_next_pixel_rgba(void* out)
    {
      colortype* pixel = (colortype*) out;
      FI_type* ptr = (FI_type*) (m_ptr);
      (*pixel)[0] = ptr[y].red;
      (*pixel)[1] = ptr[y].green;
      (*pixel)[2] = ptr[y].blue;
      (*pixel)[3] = ptr[y].alpha;

      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }

    inline
    void freeimage_reader_plugin::read_next_pixel_rgba8(void* out)
    {
      colors::rgba8* pixel = (colors::rgba8*) out;
      (*pixel)[0] = m_ptr[y * bpp + FI_RGBA_RED];
      (*pixel)[1] = m_ptr[y * bpp + FI_RGBA_GREEN];
      (*pixel)[2] = m_ptr[y * bpp + FI_RGBA_BLUE];
      (*pixel)[3] = m_ptr[y * bpp + FI_RGBA_ALPHA];
      if (++y == m_domain.pmax[1])
        {
          y = 0;
          m_ptr -= pitch;
        }
    }

    // Main method that set the right
    inline
    void freeimage_reader_plugin::read_next_line(void* out_)
    {
      mln_precondition(m_dib != NULL);
      m_read_next_line(out_);
    }

    inline
    void freeimage_reader_plugin::read_next_pixel(void* out_)
    {
      mln_precondition(m_dib != NULL);
      m_read_next_pixel(out_);
    }

    inline
    freeimage_reader_plugin::~freeimage_reader_plugin()
    {
      mln_precondition(m_dib != NULL);
      FreeImage_Unload(m_dib);
      m_dib = NULL;
    }

    /******************************************/
    /****     Writer Freeimage Plugin     ****/
    /******************************************/

    inline
    freeimage_writer_plugin::freeimage_writer_plugin(std::ostream& os,
                                                     FREE_IMAGE_FORMAT fif)
      : m_osw(os),
        m_nrows (0),
        m_ncols (0),
        m_fif(fif),
        m_dib (NULL),
        m_vtype(typeid(void))
    {
      if (m_fif == FIF_UNKNOWN || !FreeImage_FIFSupportsWriting(m_fif))
        throw MLNIOException("File format not supported for writing.");
    }



    inline
    std::function<void(void*)>
    freeimage_writer_plugin::get_write_next_line_method() const
    {
      return m_write_next_line;
    }

    inline
    std::function<void(void*)>
    freeimage_writer_plugin::get_write_next_pixel_method() const
    {
      return m_write_next_pixel;
    }

    inline
    void
    freeimage_writer_plugin::_get_fit_and_bpp_from_type(FREE_IMAGE_TYPE& fit,
                                                        unsigned& bpp,
                                                        std::type_index vtype) const
    {
      static std::unordered_map<std::type_index,
                                std::pair<FREE_IMAGE_TYPE, int> > tinfo =
        {
          { typeid(bool),   {FIT_BITMAP, 1}  },
          { typeid(uint8),  {FIT_BITMAP, 8}  },
          { typeid(rgb8),   {FIT_BITMAP, 24} },
          { typeid(colors::rgba8),  {FIT_BITMAP, 32} },
          { typeid(uint16), {FIT_UINT16, 16} },
          { typeid(int16),  {FIT_INT16,  16} },
          { typeid(uint32), {FIT_UINT32, 32} },
          { typeid(int32),  {FIT_INT32,  32} },
          { typeid(float),  {FIT_FLOAT,  32} },
          { typeid(double), {FIT_DOUBLE, 64} },
          { typeid(rgb16),   {FIT_RGB16, 48} },
          { typeid(colors::rgba16),  {FIT_RGBA16, 64} },
          { typeid(rgb<float>), {FIT_RGBF, 96} },
          { typeid(colors::rgba<float>), {FIT_RGBAF, 128} }
        };

      auto x = tinfo.find(vtype);
      if (x != tinfo.end())
        std::tie(fit, bpp) = x->second;
      else {
        fit = FIT_UNKNOWN;
        bpp = 0;
      }
    }


    inline
    bool
    freeimage_writer_plugin::can_write(std::type_index type) const
    {
      FREE_IMAGE_TYPE fit;
      unsigned bpp;
      _get_fit_and_bpp_from_type(fit, bpp, type);

      // This test fails even for legit format e.g. rgbf32
      // return FreeImage_FIFSupportsExportType(m_fif, fit) and
      //   FreeImage_FIFSupportsExportBPP(m_fif, bpp);
      return (fit != FIT_UNKNOWN);
    }


    inline
    void
    freeimage_writer_plugin::_allocate()
    {
      FREE_IMAGE_TYPE   fit;
      unsigned          red_mask = 0;
      unsigned          green_mask = 0;
      unsigned          blue_mask = 0;

      mln_precondition(m_dib == NULL);

      if (m_vtype == typeid(rgb8) or m_vtype == typeid(colors::rgba8))
        {
          red_mask = FI_RGBA_RED_MASK;
          green_mask = FI_RGBA_GREEN_MASK;
          blue_mask = FI_RGBA_BLUE_MASK;
        }

      _get_fit_and_bpp_from_type(fit, bpp, m_vtype);

      if (fit == FIT_BITMAP)
        m_dib = FreeImage_AllocateT(fit, m_ncols, m_nrows, bpp, red_mask, green_mask, blue_mask);
      else
        m_dib = FreeImage_AllocateT(fit, m_ncols, m_nrows);

      if (m_dib == NULL)
        {
          std::stringstream msg;
          msg << "Unable to allocate the image.\n"
              << "Image Type: " << fit << " BPP: " << bpp << "\n";
          throw MLNIOException(msg.str());
        }

      // Set the writing function.
      if (m_vtype == typeid(bool))
        {
          m_write_next_pixel = std::bind(&freeimage_writer_plugin::write_next_pixel_bool,
                                         this, std::placeholders::_1);
          m_write_next_line = std::bind(&freeimage_writer_plugin::write_next_line_bool,
                                         this, std::placeholders::_1);
        }
      else if (m_vtype == typeid(rgb8))
        {
          m_write_next_pixel = std::bind(&freeimage_writer_plugin::write_next_pixel_rgb8,
                                         this, std::placeholders::_1);
          m_write_next_line = std::bind(&freeimage_writer_plugin::write_next_line_rgb8,
                                         this, std::placeholders::_1);
        }
      else if (m_vtype == typeid(colors::rgba8))
        {
          m_write_next_pixel = std::bind(&freeimage_writer_plugin::write_next_pixel_rgba8,
                                         this, std::placeholders::_1);
          m_write_next_line = std::bind(&freeimage_writer_plugin::write_next_line_rgba8,
                                         this, std::placeholders::_1);
        }
      else
        {
          m_write_next_pixel = std::bind(&freeimage_writer_plugin::write_next_pixel_gray,
                                         this, std::placeholders::_1);
          m_write_next_line = std::bind(&freeimage_writer_plugin::write_next_line_gray,
                                        this, std::placeholders::_1);
        }

      // Set the palette
      RGBQUAD* palette = FreeImage_GetPalette(m_dib);
      if (m_vtype == typeid(bool))
        {
          palette[0] = {0u, 0u, 0u, 0u};
          palette[1] = {255u, 255u, 255u, 0u};
        }
      else if (m_vtype == typeid(uint8))
        {
          for (unsigned i = 0; i < 256; ++i)
            palette[i] = {(BYTE)i, (BYTE)i, (BYTE)i, 0u};
        }

      // Compute some image properties
      pitch = FreeImage_GetPitch(m_dib);
      byteperline = FreeImage_GetLine(m_dib);
      psz = byteperline / m_ncols;
      m_ptr = (char*) (FreeImage_GetBits(m_dib) + (m_nrows-1) * pitch);
      x = 0;
      y = 0;
    }


    inline
    void
    freeimage_writer_plugin::set_domain(const box2d& dom)
    {
      point2d shp = dom.shape();
      m_nrows = shp[0];
      m_ncols = shp[1];
      if (m_ncols > 0 and m_nrows > 0 and m_vtype != typeid(void))
        this->_allocate();
    }

    inline
    void
    freeimage_writer_plugin::set_value_type(std::type_index type)
    {
      m_vtype = type;
      if (m_ncols > 0 and m_nrows > 0 and m_vtype != typeid(void))
        this->_allocate();
    }

    inline
    freeimage_writer_plugin::~freeimage_writer_plugin()
    {
      FreeImageIO fio = {
        NULL,
        internal::ostream_wrapper::write,
        internal::ostream_wrapper::seek,
        internal::ostream_wrapper::tell
      };

      mln_precondition(m_dib != NULL);
      fi_handle handle = (fi_handle) (&m_osw);

      bool res = FreeImage_SaveToHandle(m_fif, m_dib, &fio, handle, 0);
      FreeImage_Unload(m_dib);
      if (!res)
        throw MLNIOException("Unable to save the image.");
    }


    /******************************************/
    /****         Writing methods         ****/
    /******************************************/

    // 1-bit DIBs are stored using each bit as an index into the color table. The most
    // significant bit is the leftmost pixel.
    inline
    void freeimage_writer_plugin::write_next_line_bool(void* src_)
    {
      mln_precondition(m_dib != NULL);

      char* src = (char*)src_;
      unsigned y, z;
      for (y = 0, z = 0; (y+8) < m_ncols; ++z)
        for (int b = 7; b >= 0; --b, ++y)
          m_ptr[z] |= (src[y] << b);

      for (int b = 7; y < m_ncols; --b, ++y)
        m_ptr[z] |= (src[y] << b);

      m_ptr -= pitch;
    }

    inline
    void freeimage_writer_plugin::write_next_pixel_bool(void* src_)
    {
      mln_precondition(m_dib != NULL);

      char* src = (char*)src_;
      m_ptr[x / 8] |= (*src) << (7-y);

      ++x;
      ++y;
      if (x == m_ncols) {
        x = 0;
        y = 0;
        m_ptr -= pitch;
      } else if (y == 8) {
        y = 0;
      }
    }

    // 8-bit DIBs are the easiest to store because each byte is an
    // index into the color table.
    // Non standard image types such as short, long, float or double do not have a color
    // table. Pixels are stored in a similar way as 8-bit DIB
    inline
    void freeimage_writer_plugin::write_next_line_gray(void* src)
    {
      std::memcpy(m_ptr, src, byteperline);
      m_ptr -= pitch;
    }

    inline
    void freeimage_writer_plugin::write_next_pixel_gray(void* src)
    {
      std::memcpy(m_ptr + y * psz, src, psz);
      if (++y == m_ncols)
        {
          y = 0;
          m_ptr -= pitch;
        }
    }


    // 24-bit DIBs have every 3 bytes representing a color, using the
    // same ordering as the RGBTRIPLE structure.
    inline
    void freeimage_writer_plugin::write_next_line_rgb8(void* src)
    {
      rgb8* buffer = (rgb8*) src;
      for (unsigned y = 0; y < m_ncols; y++)
        {
          m_ptr[y * psz + FI_RGBA_RED] =   buffer[y][0];
          m_ptr[y * psz + FI_RGBA_GREEN] = buffer[y][1];
          m_ptr[y * psz + FI_RGBA_BLUE] =  buffer[y][2];
        }
      m_ptr -= pitch;
    }

    inline
    void freeimage_writer_plugin::write_next_pixel_rgb8(void* src)
    {
      rgb8* pixel = (rgb8*) src;
      m_ptr[y * psz + FI_RGBA_RED]   = (*pixel)[0];
      m_ptr[y * psz + FI_RGBA_GREEN] = (*pixel)[1];
      m_ptr[y * psz + FI_RGBA_BLUE]  = (*pixel)[2];

      if (++y == m_ncols)
        {
          y = 0;
          m_ptr -= pitch;
        }
    }


    // 32-bit DIBs have every 4 bytes representing a color, using the
    // same ordering as the RGBQUAD structure (ordering OS-dependant)
    inline
    void freeimage_writer_plugin::write_next_line_rgba8(void* src)
    {
      colors::rgba8* buffer = (colors::rgba8*) src;
      for (unsigned y = 0; y < m_ncols; y++)
        {
          m_ptr[y * psz + FI_RGBA_RED] =   buffer[y][0];
          m_ptr[y * psz + FI_RGBA_GREEN] = buffer[y][1];
          m_ptr[y * psz + FI_RGBA_BLUE] =  buffer[y][2];
          m_ptr[y * psz + FI_RGBA_ALPHA] =  buffer[y][3];
        }
      m_ptr -= pitch;
    }

    inline
    void freeimage_writer_plugin::write_next_pixel_rgba8(void* src)
    {
      colors::rgba8* pixel = (colors::rgba8*) src;
      m_ptr[y * psz + FI_RGBA_RED]   = (*pixel)[0];
      m_ptr[y * psz + FI_RGBA_GREEN] = (*pixel)[1];
      m_ptr[y * psz + FI_RGBA_BLUE]  = (*pixel)[2];
      m_ptr[y * psz + FI_RGBA_ALPHA]  = (*pixel)[3];

      if (++y == m_ncols)
        {
          y = 0;
          m_ptr -= pitch;
        }
    }

    // Main method that set the right
    inline
    void freeimage_writer_plugin::write_next_line(void* out_)
    {
      mln_precondition(m_dib != NULL);
      m_write_next_line(out_);
    }

    inline
    void freeimage_writer_plugin::write_next_pixel(void* out_)
    {
      mln_precondition(m_dib != NULL);
      m_write_next_pixel(out_);
    }

  }

}

#endif // ! MLN_IO_FREEIMAGE_PLUGIN_HPP
