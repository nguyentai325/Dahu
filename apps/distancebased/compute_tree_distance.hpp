#ifndef COMPUTE_TREE_DISTANCE_HPP
# define COMPUTE_TREE_DISTANCE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>
# include <mln/morpho/tos/pset.hpp>
# include <mln/morpho/tos/pset_priority.hpp>
# include <mln/morpho/tos/tos.hpp>

namespace mln
{

  namespace morpho
  {

    /// \brief Compute the tree of shapes based on distance prop.
    ///
    /// Compute the tree of shapes my immerging the image in Khalimsky grid.
    /// The user is responsible to provide an image with a constant border, since the
    /// algorithm starts flooding at point ima.domain().pmin. The user may also want
    /// to perform a median interpolation to ensure correct properties on TOS.
    ///
    /// \return
    ///
    template <typename I,
	      typename Neighborhood,
	      class Distance>
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    ToSdistance(const Image<I>& ima,
		const Neighborhood& nbh,
		mln_point(I) pmin,
		Distance dist);

    /*********************************/
    /***  Implementation	  ****/
    /*********************************/
    namespace internal
    {

      template <typename V>
      void fillhole(image2d<bool>& mask, unsigned x, image2d<V>& lpar, V val)
      {
	static std::vector<unsigned> pqueue;

	pqueue.reserve(lpar.domain().size());
	pqueue.push_back(x);

	auto dindexes = wrt_delta_index(mask, c4.dpoints);
	while (!pqueue.empty())
	  {
	    x = pqueue.back();
	    pqueue.pop_back();
	    for (int n: dindexes)
	      {
		unsigned q = n + x;
		if (mask[q]) {
		  pqueue.push_back(q);
		  mask[q] = false;
		  lpar[q] = val;
		}
	      }
	  }
      }

      template <typename V>
      struct projection;



      template <>
      struct projection<rgb8>
      {
	rgb8
	operator () (const tos::irange<rgb8>& rng, rgb8 x) const
	{
	  for (int i = 0; i < 3; ++i)
	    if (x[i] < rng.lower[i]) x[i] = rng.lower[i];
	    else if (x[i] > rng.upper[i]) x[i] = rng.upper[i];
	  return x;
	}
      };

    }



    template <typename I,
	      typename Neighborhood,
	      class Distance>
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    ToSdistance(const Image<I>& ima_,
		const Neighborhood& nbh,
		mln_point(I) pmin,
		Distance dist)
    {
      using namespace mln::morpho::tos;

      typedef mln_value(I) V;
      typedef irange<V> R;
      typedef typename I::size_type size_type;

      static constexpr size_type UNPROCESSED = value_traits<size_type>::max();
      static constexpr size_type PROCESSED = 0;
      static constexpr size_type INQUEUE = 1;

      const I& ima = exact(ima_);


      // f: image of interval in Khalimsky space
      // K: image of value in Khalimsky that tells at which a level a point is inserted
      mln_ch_value(I, R) f = tos::internal::immerse(ima, productorder_less<V> ());
      mln_concrete(I) K;
      mln_ch_value(I, size_type) parent, zpar;
      mln_concrete(I) lpar;

      std::vector<size_type> S;

      S.reserve(f.domain().size());
      resize(K, f);
      resize(lpar, f);
      resize(parent, f).init(UNPROCESSED);
      extension::fill(parent, PROCESSED);

      auto dindexes = wrt_delta_index(f, nbh.dpoints);

      // Step 1: propagation
      {
	std::vector<size_type> W;
	std::vector<unsigned> vdist;
	W.reserve(f.domain().size());
	vdist.resize(f.domain().size());

	size_type p = f.index_of_point(pmin);
	W.push_back(p);
	parent[p] = PROCESSED;
	K[p] = f[p].lower;

	image2d<bool> mask;
	resize(mask, f);
	extension::fill(mask, false);

	morpho::internal::projection<V> proj;

	while (!W.empty())
	  {
	    // retrieve the closest point in the que from the current level
	    {
	      std::transform(W.begin(), W.end(), vdist.begin(),
			     [dist, proj, &lpar, &f](unsigned x) { return dist(proj(f[x], lpar[x]), lpar[x]); });
	      auto it = std::min_element(vdist.begin(), vdist.begin() + W.size());
	      int i = it - vdist.begin();
	      p = W[i];
	      W[i] = W.back();
	      W.pop_back();

	      K[p] = (vdist[i] == 0) ? lpar[p] : proj(f[p], lpar[p]);
	    }

	    V curlevel = K[p];
	    S.push_back(p);
	    copy(parent != PROCESSED, mask);
	    morpho::internal::fillhole(mask, p, lpar, curlevel);
	    parent[p] = PROCESSED;

	    mln_foreach (int k, dindexes)
	      {
		size_type q = p + k;
		if (parent[q] == UNPROCESSED)
		  {
		    parent[q] = INQUEUE;
		    //std::cout << "Insert:" << q << " @ " << K[q] << std::endl;
		    W.push_back(q);
		  }
	      }
	  }
      }

      // 2nd step: union-find

      resize(zpar, parent).init(UNPROCESSED);
      extension::fill(zpar, UNPROCESSED);

      auto is_face_2 = [](const point2d& p) { return p[0] % 2 == 0 and p[1] % 2 == 0; };

      std::equal_to<V> eq;

      int spos = S.size()-1;
      for (int i = S.size()-1; i >= 0; --i)
	{
	  size_type p = S[i];
	  parent[p] = p;
	  zpar[p] = p;

	  size_type rp = p;
	  bool face2 = is_face_2(K.point_at_index(p));

	  mln_foreach (int k, dindexes)
	    {
	      size_type q = p + k;
	      if (zpar[q] != UNPROCESSED)
		{
		  size_type r = morpho::internal::zfind_root(zpar, q);
		  if (r != rp) { // MERGE r and p
		    if (eq(K[p], K[r]) and !face2 and is_face_2(K.point_at_index(r)))
		      {
			parent[rp] = r;
			zpar[rp] = r;
			face2 = true;
			S[spos--] = rp;
			rp = r;
		      }
		    else
		      {
			parent[r] = rp;
			zpar[r] = rp;
			S[spos--] = r;
		      }
		  }
		}
	    }
	}
      S[0] = morpho::internal::zfind_root(zpar, S[0]);

      // 3rd step: canonicalization
      for (size_type p : S)
	{
	  size_type q = parent[p];
	  if (eq(K[q], K[parent[q]]))
	    parent[p] = parent[q];
	  //mln_assertion(K1::is_face_2(K.point_at_index(parent[p])));
	}


      mln_postcondition(S.size() == K.domain().size());
      mln_postcondition(S.size() == parent.domain().size());
      // All done !
      return std::make_tuple(std::move(K), std::move(parent), std::move(S));
    }

  }

}





# ifndef MLN_INCLUDE_ONLY

# endif // ! MLN_INCLUDE_ONLY

#endif // ! COMPUTE_TREE_DISTANCE_HPP
