#ifndef CUTS_HPP
# define CUTS_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/labeling/blobs.hpp>
# include <mln/accu/accumulators/infsup.hpp>
# include <mln/morpho/saturate.hpp>
# include "shape.hpp"


template<typename R, typename V,
	 class NbhFg,
	 class NbhBg,
	 typename shape_t,
	 typename Compare>
void cut_and_get_shapes(const mln::image2d<R>&				ima,
			const NbhFg&					nbhBg,
			const NbhBg&					nbhFg,
			V						lambda,
			Compare						cmp,
			shape_set<shape_t>&				shapes)
{
  using namespace mln;

  typedef short L;

  // Thresholding
  auto input = imtransform(ima, [cmp, lambda](const R& x) {
      return cmp(x, lambda);
    });


  // Labelisation
  L nlabel;
  image2d<L> imlabel;
  std::tie(imlabel, nlabel) = labeling::blobs(input, nbhFg, L());


  // Shape computation
  image2d<bool> cc;
  resize(cc, imlabel);
  accu::accumulators::infsup< point2d, productorder_less<point2d> > bbox;

  // format(std::cout, lambda) << std::endl;
  // io::imprint(input);
  // io::imprint(imlabel);

  shape_t shp(lambda, cmp, ima.domain().size(), ima.ncols());

  for (int l = 1; l <= nlabel; ++l)
    {
      morpho::saturate(imlabel == l, nbhBg, point2d(0,0), cc);

      shp.init();
      shp.set_level(lambda, cmp);
      mln_foreach(const point2d& p, where(cc))
	shp.add_point(p);

      // insert new shape in the set
      auto res = shapes.insert(std::move(shp));

      // If shape already there, update with new information
      if (not res.second)
	res.first->update_with(shp);

      mln_postcondition(res.first->islower() or res.first->isupper());
    }
}


#endif // ! CUTS_HPP
