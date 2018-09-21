#ifndef MLN_MORPHO_TOS_IMMERSE_HPP
# define MLN_MORPHO_TOS_IMMERSE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/io/imprint.hpp>

namespace mln
{

  namespace morpho
  {

    namespace tos
    {

      namespace internal
      {

	template <typename I, class Compare = std::less<mln_value(I)> >
	mln_ch_value(I, irange<mln_value(I)> )
	  immerse(const Image<I>& ima_, Compare cmp = Compare())
	{
	  static_assert(std::is_same<typename I::domain_type, box2d>::value,
			"Only box2d handled");

	  typedef mln_value(I) V;
	  typedef irange<V>    R;
	  typedef point2d      P;
	  const I& ima = exact(ima_);

	  box2d dom = ima.domain();
	  dom.pmin = dom.pmin * 2;
	  dom.pmax = dom.pmax * 2 - 1;
	  image2d<R> out(dom);

	  mln_foreach(point2d p, ima.domain())
	    {
	      V a = ima.at(p),
		b = ima.at(p + P{0,1}),
		c = ima.at(p + P{1,0}),
		d = ima.at(p + P{1,1});

	      V min1 = inf(a,b, cmp), min2 = inf(a,c, cmp);
	      V max1 = sup(a,b, cmp), max2 = sup(a,c, cmp);
	      V min3 = inf(d, inf(c, min1, cmp), cmp);
	      V max3 = sup(d, sup(c, max1, cmp), cmp);

	      point2d q = 2 * p;
	      out.at(q) = ima.at(p);
	      out.at(q + P{0,1}) = R{min1, max1};
	      out.at(q + P{1,0}) = R{min2, max2};
	      out.at(q + P{1,1}) = R{min3, max3};
	    }

	  return out;
	}

      }

    }

  }

}

#endif // ! MLN_MORPHO_TOS_IMMERSE_HPP
