#ifndef CURVATURE_HPP
# define CURVATURE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/extension/fill.hpp>

namespace mln
{


  template <typename V>
  image2d<float>
  curvature(const image2d<V>& ima)
  {
    typedef decltype(V() + V()) Vec;

    extension::fill(ima, literal::zero);

    auto norm = [] (const Vec& x) -> int { return std::abs(x[0] + x[1] + x[2]); };
    auto sqr = [] (int x) { return x*x; };

    image2d<float> curv;
    resize(curv, ima);


    mln_foreach(const point2d& p, ima.domain())
      {
	int uxx = norm(ima.at(p + point2d{0,1}) + ima.at(p + point2d{0,-1}) - 2 * ima(p));
	int uyy = norm(ima.at(p + point2d{1,0}) + ima.at(p + point2d{-1,0}) - 2 * ima(p));
	int uxy = norm(-ima.at(p + point2d{-1,-1}) + ima.at(p + point2d{-1,1}) +
			ima.at(p + point2d{1,-1}) - ima.at(p + point2d{1,1}));

	int ux = norm((ima.at(p + point2d{0,1}) - ima.at(p + point2d{0,-1})) / 2);
	int uy = norm((ima.at(p + point2d{1,0}) - ima.at(p + point2d{-1,0})) / 2);

	int den = (sqr(ux) + sqr(uy));
	if (den != 0)
	  curv(p) = std::abs((uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den)));
	else
	  curv(p) = 0;
	//curv(p) = den;
      }

    return curv;
  }

}


#endif // ! CURVATURE_HPP
