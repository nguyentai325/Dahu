#ifndef MLN_CORE_ALGORITHM_CLONE_HPP
# define MLN_CORE_ALGORITHM_CLONE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/algorithm/copy.hpp>

/// \file
/// \brief Header file for clone algorithm.

namespace mln
{

  /// \ingroup Algorithms
  /// \brief Make and return a deep copy of an image.
  ///
  /// If the image is not already concrete (e.g. a morpher with no
  /// proper memory), the image is concretized.
  ///
  /// \tparam I The input image type
  /// \param ima The image to clone.
  template <typename I>
  mln_concrete(I)
  clone(const Image<I>& ima);

  /******************/
  /* Implem         */
  /******************/

  namespace impl
  {

    // Clone by copy construction
    template <typename Image>
    mln_concrete(Image)
    clone(const Image& ima, std::true_type _use_copy_construction_)
    {
      (void) _use_copy_construction_;
      mln_concrete(Image) x(ima);
      return x;
    }

    // Clone by deep copy
    template <typename Image>
    mln_concrete(Image)
    clone(const Image& ima, std::false_type _use_copy_construction_)
    {
      (void) _use_copy_construction_;
      mln_concrete(Image) x = imconcretize(ima);
      copy(ima, x);
      return x;
    }

  }

  template <typename I>
  mln_concrete(I)
  clone(const Image<I>& ima)
  {
    return impl::clone(exact(ima), check_t<std::is_convertible<I, mln_concrete(I)>::value and
                       not image_traits<mln_concrete(I)>::shallow_copy::value> ());
  }


}

#endif // !MLN_CORE_ALGORITHM_CLONE_HPP
