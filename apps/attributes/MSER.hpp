#ifndef APPS_ATTTRIBUTES_MSER_HPP
# define APPS_ATTTRIBUTES_MSER_HPP


# include <mln/core/image/image2d.hpp>
# include <mln/core/trace.hpp>
# include <cstdlib>

enum eMSER_attribute { MSER_DIFF, MSER_RATIO, MSER_NORM };


///
///
/// \param f
/// \param K
/// \param parent
/// \param S
/// \param eps Threshold value to select the ancestor with wich the delta area will be computed
/// \param amser Attribute type:
///    *  MSER_DIFF for the delta area, \$ f(q) - f(p) \$,
///    *  MSER_RATIO for the ratio of the areas (between 0 and 1) \$ f(p) / f(q) \$
///    *  MSER_NORM for the delta area but normalized with the current area. \$ \frac{(f(q) - f(p))}{f(p)} \$
/// \param dist The distance function between two elements V (defaults to the absolute value)
/// \return
///
template <typename V, typename T, class Distance>
mln::image2d<float>
compute_MSER_attribute(const mln::image2d<V>& f,
		       const mln::image2d<T>& K,
		       const mln::image2d<unsigned>& parent,
		       const std::vector<unsigned>& S,
		       typename std::result_of<Distance(V,V)>::type eps,
		       eMSER_attribute amser = MSER_DIFF,
		       Distance dist = Distance() );

template <typename V, typename T, typename T2>
mln::image2d<float>
compute_MSER_attribute(const mln::image2d<V>& f,
		       const mln::image2d<T>& K,
		       const mln::image2d<unsigned>& parent,
		       const std::vector<unsigned>& S,
		       T2 eps,
		       eMSER_attribute amser = MSER_DIFF);



/********************/
/** Implementation **/
/********************/

template <typename V, typename T, typename T2>
mln::image2d<float>
compute_MSER_attribute(const mln::image2d<V>& f,
		       const mln::image2d<T>& K,
		       const mln::image2d<unsigned>& parent,
		       const std::vector<unsigned>& S,
		       T2 eps,
		       eMSER_attribute amser)
{
  auto fun = [](V x, V y) { return std::abs(x - y); };
  return compute_MSER_attribute(f, K, parent, S, eps, amser, fun);
}


template <typename V, typename T, class Distance>
mln::image2d<float>
compute_MSER_attribute(const mln::image2d<V>& f,
		       const mln::image2d<T>& K,
		       const mln::image2d<unsigned>& parent,
		       const std::vector<unsigned>& S,
		       typename std::result_of<Distance(V,V)>::type eps,
		       eMSER_attribute amser,
		       Distance dist)
{
  using namespace mln;
  trace::entering("compute_MSER_attribute");


  // compute area
  image2d<unsigned> area;
  image2d<unsigned> maxarea; // the max area of the children in Δh

  {
    resize(area, K).init(0);
    resize(maxarea, K).init((unsigned)0);
    area[S[0]] = 1;
    for (int i = S.size()-1; i > 0; --i)
      {
	unsigned x = S[i];
	++area[x];
	area[parent[x]] += area[x];

	if (K[x] != K[parent[x]]) { // canonical propagate to parent
	  unsigned y = parent[x];
	  unsigned a = area[x];
	  V v = f[x];
	  while (parent[y] != y and dist(v, f[y]) < eps) {
	    y = parent[y];
	  }
	  maxarea[y] = std::max(a, maxarea[y]);
	}
      }
  }

  // MSER attr
  image2d<float> mser;
  resize(mser, K);
  {
    //mser[S[0]] = 0; // FIXME !!
    for (int i = S.size()-1; i > 0; --i)
      {
	unsigned x = S[i];
	if (K[parent[x]] == K[x]) // not a canonical element
	  continue;

	V v = f[x];
	unsigned y = parent[x];
	unsigned smally = parent[x];
	while (dist(v, f[y]) < eps and y != parent[y]) {
	  if (dist(v,f[y]) < (eps / 3))
	    smally = y;
	  y = parent[y];
	}

	switch (amser) {
	  case MSER_DIFF: mser[x] = area[y] - maxarea[x]; break;
	  case MSER_RATIO: mser[x] = //(float) area[x] / area[y]; break;
	    //mser[x] = std::max(10.0f, (float) (area[y] - maxarea[x]) / area[x]); break;
	    mser[x] = (float) (area[smally] - area[x]) / (area[y] - area[x]); break;
	  case MSER_NORM: mser[x] = (float) (area[y] - maxarea[x]) / maxarea[x]; break;
	}
      }

    if (amser == MSER_DIFF or amser == MSER_NORM)
      mser[S[0]] = 0;
    else
      mser[S[0]] = 1;
  }

  trace::exiting();
  return mser;
}


#endif // ! APPS_ATTTRIBUTES_MSER_HPP
