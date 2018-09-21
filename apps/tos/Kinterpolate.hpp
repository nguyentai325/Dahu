#ifndef KINTERPOLATE_HPP
# define KINTERPOLATE_HPP

# include <mln/core/image/image2d.hpp>
# include <apps/tos/addborder.hpp>
# include <exception>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/range.hpp>



# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/fill.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>
# include <mln/morpho/tos/pset.hpp>
# include <mln/morpho/tos/pset_priority.hpp>
# include <mln/morpho/maxtree/maxtree.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/core/image/morphers/casted_image.hpp>




namespace mln
{

  /// \brief perform a median interpolation to make the iamge well-composed
  /// Note: the values are doubled to prevent rounding error
  template <typename U, typename V = decltype( std::declval<U>() * 2),
	    class Compare = std::less<U> >
  image2d<V>
  interpolate_median(const image2d<U>& ima, const V& = V(), const Compare& = Compare());


  /// \brief perform a interpolation by setting the mean on
  /// on 0-1 face from 2-faces.
  /// \param ima Original image
  /// \return An image twice as big as \p ima with mean on 0 abd 1 faces. 
  template <typename T>
  image2d<T>
  interpolate_k1(const image2d<T>& ima);

  /// \brief perform a immersion, 0-1 face have undedfined values.
  /// \param ima Original image
  /// \return An image twice as big as \p ima.
  template <typename I>
  image2d<mln_value(I)>
  immerse_k1(const Image<I>& ima, mln_value(I) v = mln_value(I) ());

  template <typename T>
  image2d<T>
  unimmerse_k1(const image2d<T>& ima);

  /// \brief Adjust ima to the domain using K1 or K2 interpolation
  /// And adding a border if necessary
  template <typename T>
  image2d<T>
  Kadjust_to(const image2d<T>& ima, box2d domain, const std::string& method = "");


  /*******************************/
  /***    Implementation       ***/
  /*******************************/

  template <typename U, typename V, class Compare>
  image2d<V>
  interpolate_median(const image2d<U>& ima, const V&, const Compare& cmp)
  {
    image2d<V> out(2*ima.nrows()-1, 2*ima.ncols()-1);

    auto med4x2 = [&cmp] (const U& a, const U& b, const U& c, const U& d) {
      U min1, min2, max1, max2;
      std::tie(min1, max1) = std::minmax(a,b, cmp);
      std::tie(min2, max2) = std::minmax(c,d, cmp);
      return std::min(max1, max2, cmp) + std::max(min1, min2, cmp);
    };

    mln_foreach(const point2d& p, ima.domain())
      {
	U a = ima.at(p),
	  b = ima.at(p + point2d{0,1}),
	  c = ima.at(p + point2d{1,0}),
	  d = ima.at(p + point2d{1,1});

	point2d q = 2 * p;
	out.at(q) = 2 * ima.at(p);
	out.at(q + point2d{0,1}) = (a + b);
	out.at(q + point2d{1,0}) = (a + c);
	out.at(q + point2d{1,1}) = med4x2(a,b,c,d);
      }
    return out;
  }


  template <typename T>
  image2d<T>
  interpolate_k1(const image2d<T>& ima)
  {
    image2d<T> out(2*ima.nrows()-1, 2*ima.ncols()-1);
    typedef point2d P;
    mln_foreach(point2d p, ima.domain())
      {
	T a = ima.at(p),
	  b = ima.at(p + P{0,1}),
	  c = ima.at(p + P{1,0}),
	  d = ima.at(p + P{1,1});

	point2d q = 2 * p;
	out.at(q) = ima.at(p);
	out.at(q + P{0,1}) = (a + b) / 2;
	out.at(q + P{1,0}) = (a + c) / 2;
	out.at(q + P{1,1}) = (a + b + c + d) / 4;
      }

    return out;
  }


  template <typename T>
  image2d<T>
  interpolate_k2(const image2d<T>& ima)
  {
    
    image2d<T> out(2*ima.nrows()-1, 2*ima.ncols()-1);
    typedef point2d P;
    mln_foreach(point2d p, ima.domain())
      {
	T a = ima.at(p),
	  b = ima.at(p + P{0,1}),
	  c = ima.at(p + P{1,0}),
	  d = ima.at(p + P{1,1});
	  
	T min1, min2, min3;
	for (int i = 0; i< 3; ++i)
	  {
		min1[i] = std::min(a[i],b[i]);
		min2[i] = std::min(a[i],c[i]);
		min3[i] = std::min(std::min(std::min(a[i],b[i]),c[i]),d[i]);		
	  }

	point2d q = 2 * p;
	out.at(q) = ima.at(p);
	out.at(q + P{0,1}) = min1;
	out.at(q + P{1,0}) = min2;
	out.at(q + P{1,1}) = min3;
      }

    return out;
  }








  template <typename I>
  image2d<mln_value(I)>
  immerse_k1(const Image<I>& ima_, mln_value(I) v)
  {
    static_assert( std::is_convertible<typename I::domain_type, box2d>::value,
		   "Image domain must be convertible to box2d." );

    const I& ima = exact(ima_);

    box2d dom = ima.domain();
    dom.pmin = dom.pmin * 2;
    dom.pmax = dom.pmax * 2 - 1;

    image2d<mln_value(I)> out(dom, 3, v);

    mln_foreach(const point2d& p, ima.domain())
      out(2*p) = ima(p);
    return out;
  }

  template <typename T>
  image2d<T>
  unimmerse_k1(const image2d<T>& ima)
  {
    sbox2d dom(ima.domain().pmin, ima.domain().pmax, point2d{2,2});
    image2d<T> out((ima.nrows() + 1) / 2, (ima.ncols() + 1) / 2);

    copy(ima | dom, out);
    return out;
  }


  template <typename T>
  image2d<T>
  Kadjust_to(const image2d<T>& ima, box2d domain, const std::string& method)
  {
    point2d shp0 = ima.domain().shape();
    point2d shp = domain.shape();
    sbox2d subdomain;

    std::function< image2d<T>(const image2d<T>&) > callback;
    if (method == "zero")
      callback = [] (const image2d<T>& ima) { return immerse_k1(ima); };
    else
      callback = &interpolate_k1<T>;


    if (shp+2 == shp0) { // remove border
      subdomain = sbox2d{ima.domain().pmin + point2d{1,1},
                         ima.domain().pmax - point2d{1,1},
                         {1,1}};
    } else if (shp*2-1 == shp0) { // unimmerse k1
      subdomain = sbox2d{ima.domain().pmin, ima.domain().pmax, {2,2}};
    } else if (shp*2+3 == shp0) { // unimmerse k1 - border
      subdomain = sbox2d{ima.domain().pmin + point2d{2,2},
                         ima.domain().pmax - point2d{2,2},
                         {2,2}};
    } else if (shp*4-3 == shp0) { // unimmerse k2
      subdomain = sbox2d{ima.domain().pmin, ima.domain().pmax, {4,4}};
    } else if (shp*4+5 == shp0) { // unimmerse k2 - border
      subdomain = sbox2d{ima.domain().pmin + point2d{4,4},
                         ima.domain().pmax - point2d{4,4},
                         {4,4}};
    } else if (shp == shp0) {
        return ima;
    } else if (shp == shp0+2) {
        return addborder(ima, lexicographicalorder_less<T>());
    } else if (shp == shp0*2-1) { // immerse_k1
      return callback(ima);
    } else if (shp == (shp0*2+3)) { // addborder + callback
      return callback(addborder(ima, lexicographicalorder_less<T>()));
    } else if (shp == (shp0*4-3)) { // immerse_k2
      return callback(callback(ima));
    } else if (shp == (shp0*4+5)) { // addborder + immerse_k2
      return callback(callback(addborder(ima, lexicographicalorder_less<T>())));
    } else {
      std::cerr << "Unable to convert the image from: " << ima.domain() << " to " << domain << "\n";
      throw std::runtime_error("Domains have invalid size.");
    }

    image2d<T> out(domain);
    copy(ima | subdomain, out);
    return out;
  }

}


#endif // ! KINTERPOLATE_HPP
