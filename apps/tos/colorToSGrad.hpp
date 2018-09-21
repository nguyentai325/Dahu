#ifndef COLORTOSGRAD_HPP
# define COLORTOSGRAD_HPP

# include <vector>
# include <mln/core/image/image2d.hpp>
# include <mln/core/colors.hpp>

namespace mln
{


  /// \brief Compute the ToS on the areas of R,G,B channel guided by gradient
  /// The algorithm internally add a border to the image.
  /// 
  /// \param ima The original ima (rgb)
  /// \param[out] K
  /// \param[out] parent
  /// \param[out] S
  ///
  /// The following step are done:
  /// * adding border
  /// * computing 3 ToS on r,g,b
  /// * computing area attributes + merge by gradient
  /// * computing ToS on area attributes.
  void colorToSGrad(const image2d<rgb8>& ima,
		    image2d<unsigned>& K,
		    image2d<unsigned>& parent,
		    std::vector<unsigned>& S);

  /// \brief Compute the ToS on the areas of R,G,B channel guided by gradient
  /// The difference with colorToSGrad is that it computes the 2nd ToS on the
  /// 2 faces only such that finally K's domain is twice as large as original domain
  /// (instead of 4X ine the case of colorToSGrad)
  /// \param ima The original ima (rgb)
  /// \param[out] K
  /// \param[out] parent
  /// \param[out] S
  void colorToSGrad_2f(const image2d<rgb8>& ima,
		       image2d<unsigned>& K,
		       image2d<unsigned>& parent,
		       std::vector<unsigned>& S);

  /// \brief Compute the ToS on the areas of R,G,B channel guided by gradient
  /// \param ima The original ima (rgb)
  /// \param[out] K
  /// \param[out] parent
  /// \param[out] S
  void colorToSGrad_with_mintree(const image2d<rgb8>& ima,
				 image2d<unsigned>& K,
				 image2d<unsigned>& parent,
				 std::vector<unsigned>& S);

}


#endif // ! COLORTOSGRAD_HPP
