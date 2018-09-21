#ifndef MLN_CORE_ALGORITHM_PASTE_HPP
# define MLN_CORE_ALGORITHM_PASTE_HPP

# include <mln/core/concept/image.hpp>

namespace mln {

  /// \brief \p paste copy each pixel from the input image to the output image.
  ///
  /// \ingroup Algorithms
  ///
  /// Each pixel of the input image is paste to the output image. It implies
  /// that \p input's domain must be included in \p output's domain.

  /// This is equivalent to the following code:
  ///
  /// \code
  /// mln_foreach (auto px, input.pixels())
  ///    output(px.point()) = px.val()
  /// \endcode
  ///
  /// \param[in] input The input image.
  /// \param[out] output The output image
  ///
  /// \return The output image.
  ///
  /// \pre `output.domain()` includes `input.domain()`
  /// \tparam InputImage is a model of the Forward Image.
  /// \tparam OutputImage is a model of the Writable Point-Accessible Image
  template <typename InputImage, typename OutputImage>
  OutputImage&
  paste(const Image<InputImage>& input, Image<OutputImage>& output);

  /// \ingroup Algorithms
  /// \overload
  template <typename InputImage, typename OutputImage>
  OutputImage&&
  paste(const Image<InputImage>& input, Image<OutputImage>&& output);


/******************************************/
/****          Implementation          ****/
/******************************************/

  namespace impl
  {
    template <typename I, typename J>
    void paste(const I& ima, J& out)
    {
      mln_pixter(px, ima);
      mln_forall(px)
        out(px->point()) = px->val();
    }

  }

  template <typename InputImage, typename OutputImage>
  OutputImage&
  paste(const Image<InputImage>& input, Image<OutputImage>& output_)

  {
    impl::paste(exact(input), exact(output_));
    return exact(output_);
  }

  template <typename InputImage, typename OutputImage>
  OutputImage&&
  paste(const Image<InputImage>& input_, Image<OutputImage>&& output_)
  {
    paste(input_, output_);
    return move_exact(output_);
  }


} // end of namespace mln

#endif //!MLN_CORE_ALGORITHM_IMFILL_HPP
