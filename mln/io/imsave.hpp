#ifndef MLN_IO_IMSAVE_HPP
# define MLN_IO_IMSAVE_HPP

# include <mln/core/trace.hpp>
# include <mln/io/saver.hpp>
# include <mln/io/freeimage_plugin.hpp>
# include <fstream>

namespace mln
{

  namespace io
  {

    template <typename I>
    void imsave(const Image<I>& ima, std::ostream& os, FREE_IMAGE_FORMAT fif = FIF_TIFF);

    template <typename I>
    void imsave(const Image<I>& ima, const std::string& path);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    template <typename I>
    void imsave(const Image<I>& ima, std::ostream& os, FREE_IMAGE_FORMAT fif)
    {
      mln_entering("mln::io::imsave");

      freeimage_writer_plugin plugin(os, fif);
      Saver2D<I> saver;
      saver.save(ima, &plugin, false);

      mln_exiting();
    }


    template <typename I>
    void imsave(const Image<I>& ima, const std::string& path)
    {
      if (path == "-")
        imsave(ima, std::cout);
      else
        {
          FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(path.c_str());
          std::ofstream os(path, std::ios::binary);

          if (fif == FIF_UNKNOWN)
            imsave(ima, os);
          else
            imsave(ima, os, fif);
          os.close();
        }
    }


  } // end of namespace mln::io
} // end of namespace mln



#endif //!MLN_IO_IMSAVE_HPP
