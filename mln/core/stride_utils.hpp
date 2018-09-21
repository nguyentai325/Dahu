#ifndef STRIDE_UTILS_HPP
# define STRIDE_UTILS_HPP

namespace mln
{

  template <size_t ndim, typename SizeType, typename CoordType, typename DifferenceType>
  void compute_delta_strides(const SizeType* strides, const CoordType* shape, DifferenceType* out);

  template <size_t ndim, typename SizeType, typename CoordType, typename DifferenceType>
  void compute_negative_delta_strides(const SizeType* strides, const CoordType* shape, DifferenceType* out);

  /********************/
  /** Implementation  */
  /********************/

  template <size_t ndim, typename SizeType, typename CoordType, typename DifferenceType>
  inline
  void
  compute_delta_strides(const SizeType* strides, const CoordType* shape, DifferenceType* out)
  {
    SizeType offset = strides[ndim-1];
    out[ndim-1] = strides[ndim-1];
    for (int i = ndim-2; i >= 0; --i)
      {
	offset *= (shape[i+1] - 1);
	mln_assertion(strides[i] - offset > 0);
	out[i] = strides[i] - offset;
      }
  }

  template <size_t ndim, typename SizeType, typename CoordType, typename DifferenceType>
  inline
  void
  compute_negative_delta_strides(const SizeType* strides, const CoordType* shape, DifferenceType* out)
  {
    SizeType offset = strides[ndim-1];
    out[ndim-1] = -strides[ndim-1];
    for (int i = ndim-2; i >= 0; --i)
      {
	offset *= (shape[i+1] - 1);
	out[i] = -(strides[i] - offset);
	mln_assertion(out[i] < 0);
      }
  }


}

#endif // ! STRIDE_UTILS_HPP
