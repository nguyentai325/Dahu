#ifndef SALIENCY_HPP
# define SALIENCY_HPP

# include <vector>
# include <mln/core/image/image2d.hpp>
# include <mln/core/algorithm/accumulate.hpp>
# include <mln/accu/accumulators/max.hpp>
# include <apps/simplification/simplify.hpp> // for compute_depth
# include <apps/tos/topology.hpp>



/**
* \brief Compute a saliency map from an saliency attribute
*
**/
template <typename T, typename V>
mln::image2d<float>
saliencymap(const mln::image2d<V>& attr,
	    const mln::image2d<T>& K,
	    const mln::image2d<unsigned>& parent,
	    const std::vector<unsigned>& S);


/**
* \brief Helper function to compute a saliency map from tagged nodes.
*
* Let say you have computed an enery on some selected nodes. Every selected node has
* a non-zero energy and non-selected node has a zero energy. This function returns
* an image where non selected node have been collapsed.
* Thus it can then be used for \see saliencymap.
*/
template <typename T, typename V>
mln::image2d<V>
collapse_zero_nodes(const mln::image2d<V>& attr,
		    const mln::image2d<T>& K,
		    const mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S);


/***********************/
/* Implementation     **/
/***********************/

namespace internal
{

}



template <typename T, typename V>
mln::image2d<float>
saliencymap(const mln::image2d<V>& attr,
	    const mln::image2d<T>& K,
	    const mln::image2d<unsigned>& parent,
	    const std::vector<unsigned>& S)
{
  static_assert(std::is_convertible<V, float>::value,
		"Attribute value type must be convertible to float");

  using namespace mln;

  mln_precondition(attr.domain() == K.domain());
  mln_precondition(parent.domain() == K.domain());

  image2d<unsigned> depth = compute_depth(K, parent, S);

  image2d<float> out;
  resize(out, K).init(0);

  //float amax = attr[S[0]];
  float amax = accumulate(attr, accu::features::max<> ());

  auto dx = wrt_delta_index(K, c8.dpoints);


  for (int i = S.size()-1; i >= 0; --i)
    {
      unsigned x = S[i];
      if (!K1::is_face_2(K.point_at_index(x)))
	continue;

      for (int i = 0; i < 8; ++i)
	{
	  unsigned q = x + dx[i];
	  if (not K.domain().has(K.point_at_index(q)))
	    continue;

	  // x and q are neighbors thus they belongs to the same branch
	  // if q <= x (using the tree order) there's nothing to do
	  //   because they're either in the same node or q is in a lower node.

	  unsigned y = K[x] == K[parent[x]] ? parent[x] : x;
	  unsigned z = K[q] == K[parent[q]] ? parent[q] : q;

	  //std::cout << y << " -- " << z << " : " << depth[y] << " " << depth[z] <<  std::endl;
	  if (depth[z] >= depth[y])
	    continue;

	  // retrieve the max saliency from [y -> z)
	  float sal = attr[y];
	  while (parent[y] != z) {
	    y = parent[y];
	    sal = std::max(attr[y], sal);
	  }
	  out[q] = std::max(out[q], sal / amax);
	}
    }
  return out;
}

template <typename T, typename V>
void
collapse_zero_nodes(const mln::image2d<V>& attr,
		    mln::image2d<T>& K,
		    mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S)
{
  using namespace mln;

  for (unsigned x: S)
    {
      unsigned q = parent[x];
      if (attr[x] == 0)
	K[x] = K[q];

      if (K[q] == K[parent[q]])
	parent[x] = parent[q];
    }
}


#endif // ! SALIENCY_HPP
