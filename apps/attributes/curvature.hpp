#ifndef APPS_ATTRIBUTE_CURVATURE_HPP
# define APPS_ATTRIBUTE_CURVATURE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/grays.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  /// \brief Add the curvature on a the 1F.
  /// The domain of the output image is twice as big as the original one.
  template <class V>
  image2d<float>
  compute_curvature(const image2d<V>& ima);

  extern template image2d<float> compute_curvature<uint8>(const image2d<uint8>&);
  extern template image2d<float> compute_curvature<uint16>(const image2d<uint16>&);

  /***************************/
  /**** Implementation     ***/
  /***************************/

  template <class V>
  image2d<float>
  compute_curvature(const image2d<V>& ima)
  {
    static_assert(std::is_integral<V>::value,
                  "V must be an integral type.");

    trace::entering("mln::compute_curvature");

    auto domain = ima.domain();
    domain.pmin *= 2;
    domain.pmax = domain.pmax * 2 - 1;
    image2d<float> curv(domain, ima.border(), 0);


    mln_foreach(const point2d& p, ima.domain())
      {
	// Right edge
	{
	  float ux = ima.at(p + point2d{0,1}) - ima.at(p);
	  float uy =
	    (ima.at(p + point2d{1,0}) - ima.at(p + point2d{-1,0}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{-1,1})) / 4.0;

	  float uxx =
	    (ima.at(p + point2d{0,-1}) - ima.at(p) -
	     ima.at(p + point2d{0,1})  + ima.at(p + point2d{0,2})) / 2.0;

	  float uxy =
	    (-ima.at(p + point2d{1,0}) + ima.at(p + point2d{-1,0}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{-1,1})) / 2.0;


	  float uyy =
	    (ima.at(p + point2d{-1,0}) + ima.at(p + point2d{-1,1}) +
	     ima.at(p + point2d{1,0}) + ima.at(p + point2d{1,1})
	     - 2 * ima.at(p) - 2 * ima.at(p + point2d{0,1})) / 2.0;


	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{0,1};
	  if (den != 0)
	    curv.at(p_) = std::abs(uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den));
	  else
	    curv.at(p_) = 0;
        }

	// Bottom edge
	{
	  float uy = ima.at(p + point2d{1,0}) - ima.at(p);

	  float ux =
	    (ima.at(p + point2d{0,1}) - ima.at(p + point2d{0,-1}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{1,-1})) / 4.0;

	  float uyy =
	    (ima.at(p + point2d{-1,0}) - ima.at(p) -
	     ima.at(p + point2d{1,0}) + ima.at(p + point2d{2,0})) / 2.0;

	  float uxx =
	    (ima.at(p + point2d{0,-1}) + ima.at(p + point2d{1,-1})
	     - 2 * ima.at(p) - 2 * ima.at(p + point2d{1,0}) +
	     ima.at(p + point2d{0,1}) + ima.at(p + point2d{1,1})) / 2.0;

	  float uxy =
	    (ima.at(p + point2d{0,-1}) - ima.at(p + point2d{0,1}) -
	     ima.at(p + point2d{1,-1}) + ima.at(p + point2d{1,1})) / 2.0;

	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{1,0};
	  if (den != 0)
	    curv.at(p_) = std::abs(uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den));
	  else
	    curv.at(p_) = 0;
	}
      }

    trace::exiting();
    mln_postcondition(all(curv >= 0));
    return curv;
  }

}

#endif // ! APPS_ATTRIBUTE_CURVATURE_HPP
