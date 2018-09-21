#ifndef SET_MEAN_ON_NODES_HPP
# define SET_MEAN_ON_NODES_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  /// \brief Set mean of values on nodes (computed from full component, node + children).
  /// \param ima The original image (same size as K)
  /// \param K The K image as outputed by the ToS algorithm
  /// \param S The S vector as outputed by the ToS algorithm
  /// \param parent The parent image as outputed by the ToS algorithm
  /// \param pred A boolean predicate B(p) -> bool that tells which pixels are accumulated for mean computation.
  template <typename V, typename W, class Pred>
  image2d<V>
  set_mean_on_node(const image2d<V>& ima,
		   const image2d<W>& K,
		   const std::vector<unsigned>& S,
		   const image2d<unsigned>& parent,
		   const Pred& pred);



  /// \brief Set mean of values on nodes (computed from flat zone, i.e node - children).
  /// \param ima The original image (same size as K)
  /// \param K The K image as outputed by the ToS algorithm
  /// \param S The S vector as outputed by the ToS algorithm
  /// \param parent The parent image as outputed by the ToS algorithm
  /// \param pred A boolean predicate B(p) -> bool that tells which pixels are accumulated for mean computation.
  template <typename V, typename W, class Pred>
  image2d<V>
  set_mean_on_node2(const image2d<V>& ima,
		    const image2d<W>& K,
		    const std::vector<unsigned>& S,
		    const image2d<unsigned>& parent,
		    const Pred& pred);


  /******************************/
  /*** Implementation         ***/
  /******************************/

  template <typename V, typename W, class Pred>
  image2d<V>
  set_mean_on_node(const image2d<V>& ima,
		   const image2d<W>& K,
		   const std::vector<unsigned>& S,
		   const image2d<unsigned>& parent,
		   const Pred& pred)
  {
    typedef typename std::conditional<std::is_scalar<V>::value,
				      accu::accumulators::mean<V>,
				      accu::accumulators::mean<V, vec3u> >::type Acc;
    trace::entering("mln::set_mean_on_node");

    assert(ima.domain() == K.domain());

    image2d<V> mean;
    image2d<Acc> accus;

    resize(accus, ima);
    resize(mean, ima);

    // Accumulate
    {
      for (int i = S.size()-1; i > 0 ; --i)
	{
	  unsigned k = S[i];
	  if (pred(ima.point_at_index(k)))
	      accus[k].take(ima[k]);
	  accus[parent[k]].take(accus[k]);
	}
      accus[S[0]].take(ima[S[0]]);
    }

    // reconstruct
    {
      mean[S[0]] = (V) accu::extractor::mean(accus[S[0]]);
      for (unsigned k : S)
	{
	  if (K[parent[k]] != K[k])
	    mean[k] = (V) accu::extractor::mean(accus[k]);
	  else
	    mean[k] = mean[parent[k]];
	}
    }
    trace::exiting();
    return mean;
  }


  /// \brief Set mean of values on nodes (computed from full component, node + children).
  /// \param ima The original image (same size as K)
  /// \param K The K image as outputed by the ToS algorithm
  /// \param S The S vector as outputed by the ToS algorithm
  /// \param parent The parent image as outputed by the ToS algorithm
  /// \param pred A boolean predicate B(p) -> bool that tells which pixels are accumulated for mean computation.
  template <typename V, typename W, class Pred>
  image2d<V>
  set_mean_on_node2(const image2d<V>& ima,
		    const image2d<W>& K,
		    const std::vector<unsigned>& S,
		    const image2d<unsigned>& parent,
		    const Pred& pred)
  {
    typedef typename std::conditional<std::is_scalar<V>::value,
				      accu::accumulators::mean<V>,
				      accu::accumulators::mean<V, vec3u> >::type Acc;
    trace::entering("mln::set_mean_on_node");

    assert(ima.domain() == K.domain());

    image2d<V> mean;
    image2d<Acc> accus;

    resize(accus, ima);
    resize(mean, ima);

    // Accumulate
    {
      for (int i = S.size()-1; i > 0 ; --i)
	{
	  unsigned k = S[i];
	  if (pred(ima.point_at_index(k)))
	      accus[k].take(ima[k]);
	  if (K[k] == K[parent[k]]) // transmit to parent if non-canoncial element.
	    accus[parent[k]].take(accus[k]);
	}
      accus[S[0]].take(ima[S[0]]);
    }

    // reconstruct
    {
      mean[S[0]] = (V) accu::extractor::mean(accus[S[0]]);
      for (unsigned k : S)
	{
	  if (K[parent[k]] != K[k])
	    mean[k] = (V) accu::extractor::mean(accus[k]);
	  else
	    mean[k] = mean[parent[k]];
	}
    }
    trace::exiting();
    return mean;
  }


} // end of namespce mln


#endif // ! SET_MEAN_ON_NODES_HPP
