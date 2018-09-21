#ifndef MLN_CORE_WRT_OFFSET_HPP
# define MLN_CORE_WRT_OFFSET_HPP

//# include <mln/core/range.hpp>
//#include <mln/core/std/array.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/range/range.hpp>
# include <array>

namespace mln {

  template <typename I, typename SiteSet, typename STLOutputIterator>
  inline
  void wrt_offset(const Image<I>& ima_,
		  const SiteSet& dpoints,
		  STLOutputIterator out)
  {
    const I& ima = exact(ima_);
    const size_t* strides = ima.strides();
    auto it = rng::iter(dpoints);
    for (it.init(); !it.finished(); it.next(), ++out)
      {
        *out = 0;
        for (int i = 0; i < I::ndim; ++i)
          *out += strides[i] * (*it)[i];
      }
  }


  template <typename I, size_t N>
  inline
  void wrt_offset(const Image<I>& ima_,
                  const std::array<mln_point(I), N>& dpoints,
                  std::array<typename I::difference_type, N>& out)
  {
    const I& ima = exact(ima_);
    const size_t* strides = ima.strides();
    for (unsigned j = 0; j < N; ++j)
      {
        out[j] = 0;
        for (int i = 0; i < I::ndim; ++i)
          out[j] += strides[i] * dpoints[j][i];
      }
  }

  template <typename I, size_t N>
  inline
  void wrt_delta_index(const Image<I>& ima_,
		       const std::array<mln_point(I), N>& dpoints,
		       std::array<typename I::difference_type, N>& out)
  {
    const I& ima = exact(ima_);
    const size_t* strides = ima.index_strides();
    for (unsigned j = 0; j < N; ++j)
      {
        out[j] = 0;
        for (int i = 0; i < I::ndim; ++i)
          out[j] += strides[i] * dpoints[j][i];
      }
  }

  template <typename I, typename SiteSet, typename OutputIterator>
  inline
  void wrt_delta_index(const Image<I>& ima_,
		       const SiteSet& dpoints,
		       OutputIterator out)
  {
    const I& ima = exact(ima_);
    const size_t* strides = ima.index_strides();
    auto it = rng::iter(dpoints);
    for (it.init(); !it.finished(); it.next(), ++out)
      {
        *out = 0;
        for (int i = 0; i < I::ndim; ++i)
          *out += strides[i] * (*it)[i];
      }
  }


  template <typename Image, size_t N>
  std::array<typename Image::difference_type, N>
  wrt_offset(const Image& ima,
	     const std::array<typename Image::point_type, N>& dpoints)
  {
    std::array<typename Image::difference_type, N> out;
    wrt_offset(ima, dpoints, out);
    return out;
  }

  template <typename Image, size_t N>
  std::array<typename Image::difference_type, N>
  wrt_delta_index(const Image& ima,
		  const std::array<typename Image::point_type, N>& dpoints)
  {
    std::array<typename Image::difference_type, N> out;
    wrt_delta_index(ima, dpoints, out);
    return out;
  }

} // end of namespace mln


#endif //!MLN_CORE_WRT_OFFSET_HPP
