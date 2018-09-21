#ifndef THICKEN_HPP
# define THICKEN_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/io/imprint.hpp>
# include <mln/morpho/maxtree_ufind_parallel.hpp>

namespace mln
{


  template <typename V>
  image2d<unsigned>
  thicken_bup(const image2d<V>& ima, image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    image2d<unsigned> depth;
    resize(depth, ima);

    auto offset = wrt_delta_index(ima, c4.dpoints);

    depth[S[0]] = 0;
    for (unsigned p : S)
      {
	unsigned q = parent[p];

	if (ima[p] == ima[q])
	  depth[p] = depth[q];
	else
	  depth[p] = depth[q] + 1;
      }


    image2d<unsigned> mindepth; // The depth of the highest node on the contour
    image2d<unsigned> validdepth; // The maximal depth required to be valid
    resize(mindepth, ima, ima.border(), value_traits<unsigned>::max());
    resize(validdepth, ima, ima.border(), value_traits<unsigned>::max());
    for (int i = S.size()-1; i >= 0; --i)
      {
        unsigned p = S[i];
        unsigned q = ima[p] == ima[parent[p]] ? parent[p] : p;

        mln_foreach(int k, offset)
          mindepth[q] = std::min(mindepth[q], depth[p + k]);

        if (p == q) // p is a node
          {
	    q = parent[p];
            mindepth[q] = std::min(mindepth[q], mindepth[p]);
            if (depth[p] <= validdepth[p]) // Then p is a kept
              {
                validdepth[p] = depth[p];
                validdepth[q] = std::min(validdepth[q], mindepth[p]);
              }
            else // p collapsed with its parent
              validdepth[q] = std::min(validdepth[q], validdepth[p]);
          }
      }
    io::imsave(transform(depth, [](unsigned x) -> float { return x; }), "depth.tiff");
    io::imsave(transform(validdepth, [](unsigned x) -> float { return x; }), "validdepth.tiff");

    // propagate and collapse parent
    for (unsigned p: S)
      {
        unsigned q = parent[p];
        if (ima[p] == ima[q])
          validdepth[p] = validdepth[q];
	else if (depth[p] > validdepth[p])
	  validdepth[p] = validdepth[q];

        if (validdepth[q] == validdepth[parent[q]])
          parent[p] = parent[q];
      }

    //io::imprint(depth);
    //io::imprint(validdepth);
    return validdepth;
  }

  template <typename V, typename Nbh>
  image2d<unsigned>
  thicken_tdn(const image2d<V>& ima, image2d<unsigned>& parent, const std::vector<unsigned>& S, const Nbh& nbh = c4)
  {
    image2d<unsigned> depth;
    static constexpr unsigned UNDEF = value_traits<unsigned>::max();
    resize(depth, ima).init(UNDEF);

    auto offset = wrt_delta_index(ima, nbh.dpoints);
    depth[S[0]] = 0;
    for (unsigned p : S)
      {
	unsigned q = parent[p];
	if (ima[p] == ima[q])
	  depth[p] = depth[q];
	else
	  depth[p] = depth[q] + 1;
      }



    image2d<unsigned> outdepth; // The depth of the highest node in the internal hole
    resize(outdepth, ima).init(UNDEF);
    outdepth[S[0]] = 0;
    for (unsigned p : S)
      {
        unsigned q = ima[p] == ima[parent[p]] ? parent[p] : p;

	outdepth[p] = outdepth[q];

        mln_foreach(int k, offset)
	  {
	    unsigned n = p + k;

	    if (depth[n] == UNDEF) // outside domain
	      continue;

	    unsigned r = ima[parent[n]] == ima[n] ? parent[n] : n;
	    if (outdepth[r] == UNDEF)
	      outdepth[r] = outdepth[q]+1;
	    else if (depth[q] < depth[r] and outdepth[q] > outdepth[r]) // node collapse
	      {
		unsigned x = q;
		unsigned d = outdepth[r];
		while (outdepth[x] > d) {
		  outdepth[x] = d;
		  x = parent[x];
		}
	      }

	  }
      }

    for (unsigned p : S)
      {
	unsigned q = parent[p];
	if (ima[p] == ima[q])
	  outdepth[p] = outdepth[q];

	if (outdepth[q] == outdepth[parent[q]])
	  parent[p] = parent[q];
      }


    return outdepth;
  }


}

#endif // ! THICKEN_HPP
