#ifndef ADDBORDER_HPP
# define ADDBORDER_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/neighb2d.hpp>
# include <vector>

namespace mln
{

  /**
  *
  * Important note:
  * this method assume that for two values a,b s.t a < b
  * then  a < (a+b)/2 < b
  */
  template < class V, class Compare = std::less<V> >
  image2d<V>
  addborder(const image2d<V>& ima, const Compare& cmp = Compare ())
  {
    //const I& ima = exact(ima_);
    image2d<V> out(ima.nrows() + 2, ima.ncols() + 2);

    {
      box2d box = ima.domain();
      box.pmin += 1; box.pmax += 1;
      copy(ima, out | box);
    }

    V median;
    unsigned ncols = ima.ncols(), nrows = ima.nrows();
    {
      std::vector<V> border;
      border.reserve(2 * (nrows + ncols) - 4);
      for (unsigned i = 0; i < ncols; ++i) {
	border.push_back(ima.at(0,i));
	border.push_back(ima.at(nrows-1,i));
      }

      for (unsigned i = 1; i < nrows-1; ++i) {
	border.push_back(ima.at(i,0));
	border.push_back(ima.at(i,ncols-1));
      }

      std::partial_sort(border.begin(), border.begin() + border.size()/2+1, border.end(), cmp);
      if (border.size() % 2 == 0) {
	//V a = border[border.size()/2 - 1], b = border[border.size()/2];
	//median = a + (b-a) / 2;
	median = border[border.size()/2];
      } else
	median = border[border.size()/2];
    }

    {
      for (unsigned i = 0; i < ncols+2; ++i) {
	out.at(0,i) = median;
	out.at(nrows+1,i) = median;
      }

      for (unsigned i = 1; i < nrows+1; ++i) {
	out.at(i,0) = median;
	out.at(i,ncols+1) = median;
      }
    }
    return out;
  }

  // Add a border with the median computed marginally on each channel.
  template < class V>
  image2d<V>
  addborder_marginal(const image2d<V>& ima)
  {
    //const I& ima = exact(ima_);
    image2d<V> out(ima.nrows() + 2, ima.ncols() + 2);

    for (unsigned i = 0; i < value_traits<V>::ndim; ++i)
      copy(addborder(eval(channel(ima,i))), channel(out,i));

    return out;
  }


  template <class V, class M, class Compare = std::less<V> >
  std::pair< image2d<V>, image2d<bool> >
  addborder2(const image2d<V>& ima, const Image<M>& mask_, const Compare& cmp = Compare ())
  {

    const M& mask = exact(mask_);
    image2d<V> out(ima.nrows() + 2, ima.ncols() + 2);
    image2d<bool> omask(ima.nrows() + 2, ima.ncols() + 2);

    std::vector<V> border;



    mln_foreach(point2d p, ima.domain())
      if (mask(p))
        {
          mln_foreach(point2d n, c4(p)) {
            if (!mask.domain().has(n) or !mask(n))
              border.push_back(ima(p));
            omask(n + point2d{1,1}) = true;
          }
          point2d q = p + point2d{1,1};
          omask(q) = true;
          out(q) = ima(p);
        }

    std::partial_sort(border.begin(), border.begin() + border.size()/2+1, border.end(), cmp);
    V median = border[border.size()/2];

    mln_foreach(point2d p, ima.domain())
      if (mask(p))
        {
          mln_foreach(point2d n, c4(p))
            if (!mask.domain().has(n) or !mask(n))
              {
                point2d q = n + point2d{1,1};
                out(q) = median;
              }
        }

    return std::make_pair(out, omask) ;
  }

}



#endif // ! ADDBORDER_HPP
