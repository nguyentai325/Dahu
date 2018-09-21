#ifndef APPS_TOS_ROUTINES_HPP
# define APPS_TOS_ROUTINES_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/accu/accumulator.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/core/trace.hpp>
# include <vector>

/// \file
/// \brief Common routines used with ToS.


/// \brief Compute an attribute filter inplace.
///
/// K and parent image are modified to reflect the level lines
/// removal.
template <typename T, typename V>
void
attributefilter_inplace(const mln::image2d<V>& A,
			mln::image2d<T>& K,
			mln::image2d<unsigned>& parent,
			const std::vector<unsigned>& S,
			float threshold);

/// \brief Compute a grain filter inplace.
/// \param pointfilter A function P -> bool used to filter
///                    points to counted (e.g K1::is_face_2)
template <typename T, class F>
void
grainfilter_inplace(mln::image2d<T>& K,
		    mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S,
		    float threshold,
		    F pointfilter);

/// \brief Compute a grain filter inplace.
template <typename T>
void
grainfilter_inplace(mln::image2d<T>& K,
		    mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S,
		    float threshold);

/// \brief Compute attribute per node
///
/// The attribute is computed for each node but does not
/// consider the subtree of each node. The accumulation
/// is made on \p f image based on the tree (K, parent, S).
/// Thus, it makes sense for gray level tree to use f = K.
///
/// \param pointfilter Only accumulate on some points (e.g. K1::is_face_2)
/// \param eq Equivalence relation on K s.t. y = parent(x) ^ eq(K(x), K(y)) => x,y are in the same node
template <typename V, typename T, class AccuLike, class Filter,
	  class Equiv = std::equal_to<T> >
mln::image2d< typename mln::accu::result_of<AccuLike, V>::type >
compute_attribute_per_node(const mln::image2d<V>& f,
			   const mln::image2d<T>& K,
			   const mln::image2d<unsigned>& parent,
			   const std::vector<unsigned>& S,
			   const mln::AccumulatorLike<AccuLike>& accu,
			   Filter pointfilter,
			   Equiv eq = Equiv());

template <typename V, typename T, class Filter,
	  class Equiv = std::equal_to<T> >
mln::image2d< typename mln::accu::result_of<mln::accu::features::mean<>, V>::type >
compute_mean_on_node(const mln::image2d<V>& f,
		     const mln::image2d<T>& K,
		     const mln::image2d<unsigned>& parent,
		     const std::vector<unsigned>& S,
		     Filter pointfilter,
		     Equiv eq = Equiv());

template <typename V, typename T>
mln::image2d< typename mln::accu::result_of<mln::accu::features::mean<>, V>::type >
compute_mean_on_node(const mln::image2d<V>& f,
		     const mln::image2d<T>& K,
		     const mln::image2d<unsigned>& parent,
		     const std::vector<unsigned>& S);

template <typename V>
mln::image2d<unsigned>
compute_depth(const mln::image2d<V>& K,
	      const mln::image2d<unsigned>& parent,
	      const std::vector<unsigned>& S);

/// \brief Compute attribute per node
/// The attribute is computed for each node but does
/// consider the subtree of each node.
template <typename V, typename T, class AccuLike, class Filter,
	  class Equiv = std::equal_to<T> >
mln::image2d< typename mln::accu::result_of<AccuLike, V>::type >
attribute_compute(const mln::image2d<V>& f,
		  const mln::image2d<T>& K,
		  const mln::image2d<unsigned>& parent,
		  const std::vector<unsigned>& S,
		  const mln::AccumulatorLike<AccuLike>& accu,
		  Filter pointfilter,
		  Equiv eq = Equiv());



/// \brief Low Level tos manipulation routines
/// \{

template <typename V>
std::vector<unsigned>
get_real_nodes(const mln::image2d<V>& K,
	       const mln::image2d<unsigned>& parent,
	       const std::vector<unsigned>& S);



template <typename T>
unsigned
get_canonical(const mln::image2d<T>& K,
	      const mln::image2d<unsigned>& parent,
	      unsigned x);

/// \brief Least common ancestor
/// Note that depth of root SHALL be 0.
/// -1 is the neutral element for lca i.e.
/// lca(x, -1) = x
/// lca(-1, x) = x
/// \param x, y must be canonical elements
unsigned
lca(const mln::image2d<unsigned>& depth,
    const mln::image2d<unsigned>& parent,
    unsigned x, unsigned y);
/// \}

/*********************************/
/** Implementation              **/
/*********************************/


template <typename T, typename V>
void
attributefilter_inplace(const mln::image2d<V>& A,
			mln::image2d<T>& K,
			mln::image2d<unsigned>& parent,
			const std::vector<unsigned>& S,
			float threshold)
{
  using namespace mln;
  trace::entering("Attribute filter");

  for (unsigned x: S)
    {
      unsigned q = parent[x];
      unsigned r = parent[q];

      if (A[x] < threshold) {
	K[x] = K[q];
      }
      if (K[q] == K[r]) {
	parent[x] = r;
      }
    }

  trace::exiting();
}


template <typename T, class F>
void
grainfilter_inplace(mln::image2d<T>& K,
		    mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S,
		    float threshold,
		    F pointfilter)
{
  using namespace mln;

  // compute area
  image2d<unsigned> A;
  resize(A, K).init(0);
  A[S[0]] = 1;
  for (int i = S.size()-1; i > 0; --i)
    {
      unsigned x = S[i];
      if (pointfilter(K.point_at_index(x)))
	A[x] += 1;
      A[parent[x]] += A[x];
    }

  // filter
  attributefilter_inplace(A, K, parent, S, threshold);
}



template <typename T>
void
grainfilter_inplace(mln::image2d<T>& K,
		    mln::image2d<unsigned>& parent,
		    const std::vector<unsigned>& S,
		    float threshold)
{
  using namespace mln;

  // compute area
  image2d<unsigned> A;
  resize(A, K).init(1);

  for (int i = S.size()-1; i > 0; --i)
    {
      unsigned x = S[i];
      A[parent[x]] += A[x];
    }

  // filter
  attributefilter_inplace(A, K, parent, S, threshold);
}


template <typename V, typename T, class AccuLike, class Filter, class Equiv>
mln::image2d< typename mln::accu::result_of<AccuLike, V>::type >
compute_attribute_per_node(const mln::image2d<V>& f,
			   const mln::image2d<T>& K,
			   const mln::image2d<unsigned>& parent,
			   const std::vector<unsigned>& S,
			   const mln::AccumulatorLike<AccuLike>& accu_,
			   Filter pointfilter,
			   Equiv eq)
{
  using namespace mln;
  trace::entering("attribute_compute_per_node");

  typedef typename mln::accu::result_of<AccuLike, V>::type R;
  auto accu = accu::make_accumulator<AccuLike, V>(exact(accu_));

  accu.init();

  image2d<decltype(accu)> A;
  resize(A, K).init(accu);

  mln_precondition(are_indexes_compatible(f, K));
  mln_precondition(are_indexes_compatible(f, parent));
  mln_precondition(are_indexes_compatible(f, A));

  // Compute attribute
  for (int i = S.size()-1; i >= 0; --i)
    {
      unsigned x = S[i];
      unsigned q = parent[x];
      if (pointfilter(K.point_at_index(x)))
	{
	  if (eq(K[x], K[q]))
	    A[q].take(f[x]);
	  else
	    A[x].take(f[x]);
	}
    }

  // Propagate
  image2d<R> out;
  resize(out, K);
  out[S[0]] = A[S[0]].to_result();
  for (unsigned x: S)
    {
      unsigned q = parent[x];
      if (eq(K[x], K[q]))
	out[x] = out[q];
      else
	out[x] = A[x].to_result();
    }

  trace::exiting();
  return out;
}

template <typename V, typename T, class AccuLike, class Filter, class Equiv>
mln::image2d< typename mln::accu::result_of<AccuLike, V>::type >
attribute_compute(const mln::image2d<V>& f,
		  const mln::image2d<T>& K,
		  const mln::image2d<unsigned>& parent,
		  const std::vector<unsigned>& S,
		  const mln::AccumulatorLike<AccuLike>& accu_,
		  Filter pointfilter,
		  Equiv eq)
{
  using namespace mln;
  trace::entering("attribute_compute");

  typedef typename mln::accu::result_of<AccuLike, V>::type R;
  auto accu = accu::make_accumulator<AccuLike, V>(exact(accu_));

  accu.init();

  image2d<decltype(accu)> A;
  resize(A, K).init(accu);

  mln_precondition(are_indexes_compatible(f, K));
  mln_precondition(are_indexes_compatible(f, parent));
  mln_precondition(are_indexes_compatible(f, A));

  // Compute attribute
  for (int i = S.size()-1; i >= 0; --i)
    {
      unsigned x = S[i];
      unsigned q = parent[x];
      if (pointfilter(K.point_at_index(x)))
	A[x].take(f[x]);
      A[q].take(A[x]);
    }

  // Propagate
  image2d<R> out;
  resize(out, K);
  out[S[0]] = A[S[0]].to_result();
  for (unsigned x: S)
    {
      unsigned q = parent[x];
      if (eq(K[x], K[q]))
	out[x] = out[q];
      else
	out[x] = A[x].to_result();
    }

  trace::exiting();
  return out;

}




template <typename V, typename T, class Filter, class Equiv>
mln::image2d< typename mln::accu::result_of<mln::accu::features::mean<>, V>::type >
compute_mean_on_node(const mln::image2d<V>& f,
		     const mln::image2d<T>& K,
		     const mln::image2d<unsigned>& parent,
		     const std::vector<unsigned>& S,
		     Filter pointfilter,
		     Equiv eq)
{
  return compute_attribute_per_node(f, K, parent, S,
				    mln::accu::features::mean<> (),
				    pointfilter, eq);
}


template <typename V, typename T>
mln::image2d< typename mln::accu::result_of<mln::accu::features::mean<>, V>::type >
compute_mean_on_node(const mln::image2d<V>& f,
		     const mln::image2d<T>& K,
		     const mln::image2d<unsigned>& parent,
		     const std::vector<unsigned>& S)
{
  return compute_attribute_per_node(f, K, parent, S,
				    mln::accu::features::mean<> (),
				    [] (mln::point2d ) { return true; });
}



template <typename T>
mln::image2d<unsigned>
compute_depth(const mln::image2d<T>& K,
	      const mln::image2d<unsigned>& parent,
	      const std::vector<unsigned>& S)
{
  mln::image2d<unsigned> depth;
  mln::resize(depth, K);

  // Compute depth attribute
  {
    depth[S[0]] = 0;
    for (unsigned i = 1; i < S.size(); ++i)
      {
	unsigned x = S[i];
	if (K[x] != K[parent[x]]) // canonical element
	  depth[x] = depth[parent[x]] + 1;
	else
	  depth[x] = depth[parent[x]];
      }
  }
  return depth;
}



template <typename T>
inline
unsigned
get_canonical(const mln::image2d<T>& K,
	      const mln::image2d<unsigned>& parent,
	      unsigned x)
{
  unsigned q = parent[x];
  return K[x] == K[q] ? q : x;
}


template <typename V>
std::vector<unsigned>
get_real_nodes(const mln::image2d<V>& K,
	       const mln::image2d<unsigned>& parent,
	       const std::vector<unsigned>& S)
{
  std::vector<unsigned> R;
  R.reserve(S.size());

  R.push_back(S[0]);
  for (unsigned x: S)
    if (K[x] != K[parent[x]])
      R.push_back(x);


  return R;
}



#endif // !APPS_TOS_ROUTINES_HPP
