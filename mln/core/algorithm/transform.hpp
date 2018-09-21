#ifndef MLN_CORE_ALGORITHM_TRANSFORM_HPP
# define MLN_CORE_ALGORITHM_TRANSFORM_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
/// \file

namespace mln {

  /// \ingroup Algorithms
  /// \brief Transform the value of an image through a function.
  ///
  /// This is equivalent to the following code.
  /// \code
  /// mln_iter(vout, out.values())
  /// mln_iter(vin, in.values())
  /// mln_forall(vin, vout)
  ///    *vout = f(*vin)
  ///
  /// \endcode
  ///
  /// \tparam InputImage A model of :concept:ref:`forward_image`
  /// \tparam OutputImage A model of Writable :concept:ref:`forward_image`
  /// \tparam UnaryFunction A mondel of unary function.
  /// \param input The input image.
  /// \param f The unary function.
  /// \param output The output image.
  /// \see mln::imtransform
  template <typename InputImage, typename OutputImage, class UnaryFunction>
  OutputImage&
  transform(const Image<InputImage>& input, UnaryFunction f, Image<OutputImage>& output);

  /// \overload
  template <typename InputImage, typename OutputImage, class UnaryFunction>
  OutputImage&&
  transform(const Image<InputImage>& input, UnaryFunction f, Image<OutputImage>&& output);

  template <typename InputImage, class UnaryFunction>
  unary_image_expr<UnaryFunction, InputImage>
  lazy_transform(InputImage&& input, UnaryFunction f);

  /// \overload
  template <typename InputImage, class UnaryFunction>
  mln_ch_value(InputImage,
	       typename std::decay<
		 typename std::result_of<UnaryFunction(mln_value(InputImage))>::type
		 >::type)
  transform(const Image<InputImage>& input, UnaryFunction f);


/******************************************/
/****          Implementation          ****/
/******************************************/

  namespace impl
  {
    template <typename I, typename J, class UnaryFunction>
    void
    transform(const I& input, UnaryFunction f, J& output)
    {
      mln_viter(vin, vout, input, output);
      mln_forall(vin, vout)
        *vout = f(*vin);
    }

  }


  template <typename InputImage, typename OutputImage, class UnaryFunction>
  OutputImage&
  transform(const Image<InputImage>& input, UnaryFunction f, Image<OutputImage>& output)
  {
    OutputImage& out = exact(output);
    mln_entering("mln::transform");
    impl::transform(exact(input), f, exact(output));
    mln_exiting();
    return out;
  }

  template <typename InputImage, typename OutputImage, class UnaryFunction>
  OutputImage&&
  transform(const Image<InputImage>& input, UnaryFunction f, Image<OutputImage>&& output)
  {
    mln::transform(input, f, output);
    return move_exact(output);
  }


  template <typename InputImage, class UnaryFunction>
  mln_ch_value(InputImage,
	       typename std::decay<
		 typename std::result_of<UnaryFunction(mln_value(InputImage))>::type
		 >::type)
    transform(const Image<InputImage>& input, UnaryFunction f)
  {
    typedef typename std::decay<
      typename std::result_of<UnaryFunction(mln_value(InputImage))>::type
      >::type T;

    mln_ch_value(InputImage, T) out = imchvalue<T>(exact(input));
    mln::transform(input, f, out);
    return out;
  }

  template <typename InputImage, class UnaryFunction>
  unary_image_expr<UnaryFunction, InputImage>
  lazy_transform(InputImage&& input, UnaryFunction f)
  {
    return make_unary_image_expr(std::forward<InputImage>(input), f);
  }


}



#endif // ! MLN_CORE_ALGORITHM_TRANSFORM_HPP
