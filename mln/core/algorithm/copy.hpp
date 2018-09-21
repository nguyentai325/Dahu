#ifndef MLN_CORE_ALGORITHM_COPY_HPP
# define MLN_CORE_ALGORITHM_COPY_HPP

# include <mln/core/assert.hpp>
# include <mln/core/concept/image.hpp>
/// \file

namespace mln {

  /*
  // \brief Copy an image to memory pointed by the iterator \p it.
  //
  // \param[in] input Input image
  // \param it Output iterator
  //
  // \return it + number of elements in the image.
  //
  template <typename I, typename OutputIterator>
  OutputIterator
  copy(const Image<I>& input, OutputIterator it);
  */

  /// \brief Copy a source image to an output image regardless domain.
  ///
  /// \ingroup Algorithms
  ///
  /// The output image must have a domain larger than the input image.
  /// This is equivalent to the following code.
  ///
  /// \code
  /// mln_iter(vin, input.values())
  /// mln_iter(vout, output.values())
  /// mln_forall(vin, vout)
  ///   *vout = *vin;
  /// \endcode
  ///
  /// \param[in] input Input Image
  /// \param[out] output Output Image
  /// \return The image where values have been copied in.
  ///
  /// \todo add specialization for raw images
  ///
  template <typename InputImage, typename OutputImage>
  OutputImage&
  copy(const Image<InputImage>& input, Image<OutputImage>& output);


  /// \overload
  /// \ingroup Algorithms
  template <typename InputImage, typename OutputImage>
  OutputImage&&
  copy(const Image<InputImage>& input, Image<OutputImage>&& output);



  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace impl
  {

    template <typename I, typename J>
    inline
    void
    copy(const I& input, J&& output)
    {
      mln_viter(vin, vout, input, output);

      mln_forall(vin, vout)
        *vout = (mln_value(I)) *vin;
    }

  }


  template <typename InputImage, typename OutputImage>
  OutputImage&
  copy(const Image<InputImage>& input, Image<OutputImage>& output)
  {
    static_assert(std::is_convertible<mln_value(InputImage), mln_value(OutputImage)>::value,
                  "The input image value type must be convertible to the output image value type");

    impl::copy(exact(input), exact(output));
    return exact(output);
  }

  template <typename InputImage, typename OutputImage>
  OutputImage&&
  copy(const Image<InputImage>& input, Image<OutputImage>&& output)
  {
    static_assert(std::is_convertible<mln_value(InputImage), mln_value(OutputImage)>::value,
                  "The input image value type must be convertible to the output image value type");

    impl::copy(exact(input), exact(output));
    return move_exact(output);
  }


} // end of namespace mln


#endif //!MLN_CORE_ALGORITHM_COPY_HPP
