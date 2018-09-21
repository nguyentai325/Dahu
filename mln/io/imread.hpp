#ifndef MLN_IO_IMREAD_HH
# define MLN_IO_IMREAD_HH

# include <mln/core/trace.hpp>
# include <mln/io/loader.hpp>
# include <mln/io/freeimage_plugin.hpp>
# include <fstream>

namespace mln
{
  namespace io
  {

    /// \brief Load an image from disk.
    /// \ingroup io
    ///
    /// \param path Path to the image
    /// \param[out] out Image to be loaded
    template <typename I>
    void imread(const std::string& path, Image<I>& out, bool permissive = false);

    template <typename I>
    void imread(std::istream& path, Image<I>& out, bool permissive = false);


    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    template <typename I>
    void imread(std::istream& stream, Image<I>& out, bool permissive)
    {
      mln_entering("mln::io::imread");
      freeimage_reader_plugin plugin;
      Loader2D<I> loader;
      loader.load(stream, out, &plugin, permissive);
      mln_exiting();
    }


    template <typename I>
    void imread(const std::string& path, Image<I>& out, bool permissive)
    {
      mln_entering("mln::io::imread");
      freeimage_reader_plugin plugin;
      Loader2D<I> loader;
      loader.load(path, out, &plugin, permissive);
      mln_exiting();
    }


  } /*end of namespace mln::io */
} /* end of namespace mln */
# endif
