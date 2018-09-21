#ifndef MAXTREE_UFIND_RANK_HPP
# define MAXTREE_UFIND_RANK_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/sort_indexes.hpp>
# include <mln/core/wrt_offset.hpp>
//# include <tbb/parallel_reduce.h>
# include <mln/io/imprint.hpp>
namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      static
      inline
      std::size_t
      zfindroot(image2d<std::size_t>& parent, std::size_t p)
      {
        if (parent[p] == p)
          return p;
        else
          return parent[p] = zfindroot(parent, parent[p]);
      }

    }

    template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
    std::pair< image2d<std::size_t>, std::vector<std::size_t> >
    maxtree_ufindbyrank(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
    {
      image2d<std::size_t> parent, zpar, root;
      image2d<unsigned> rank;
    image2d<bool> deja_vu;
    resize(parent, ima);
    resize(zpar, ima);
    resize(root, ima);
    resize(rank, ima, ima.border(), 0);
    resize(deja_vu, ima, ima.border(), false);

    extension::fill(deja_vu, false);

    std::vector<std::size_t> S = sort_indexes(ima, cmp);
    auto offsets = wrt_delta_index(ima, nbh.dpoints);

    for (int i = S.size()-1; i >= 0; --i)
      {
	std::size_t p = S[i];
        //std::cout << "Processing:" << p << " @ " << (int) ima[p] << std::endl;
	// make set
	{
	  parent[p] = p;
	  zpar[p] = p;
          root[p] = p;
	  deja_vu[p] = true;
	}

        std::size_t x = p; // zpar of p
	for (unsigned k = 0; k < offsets.size(); ++k)
	  {
	    std::size_t q = p + offsets[k];
	    if (deja_vu[q])
	      {
                std::size_t r = internal::zfindroot(zpar, q);
		if (r != x) // make union
		  {
                    //std::cout << "Mergin: " << root[r] << "->" << p << std::endl;
                    parent[root[r]] = p;
                    if (rank[x] < rank[r]) { //we merge p to r
                      zpar[x] = r;
                      root[r] = p;
                      x = r;
                    } else if (rank[r] < rank[p]) { // merge r to p
                      zpar[r] = p;
                    } else { // same height
                      zpar[r] = p;
                      rank[p] += 1;
                    }
		  }
	      }
	  }
        //io::imprint(parent);
      }

    // canonization
    for (std::size_t p: S)
      {
	std::size_t q = parent[p];
	if (ima[parent[q]] == ima[q])
	  parent[p] = parent[q];
      }

    return std::make_pair(std::move(parent), std::move(S));
  }

  }

}

#endif // !MLN_MORPHO_MAXTREE_UFIND_RANK_HPP
