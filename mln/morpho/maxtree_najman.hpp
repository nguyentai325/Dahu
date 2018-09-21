#ifndef MAXTREE_NAJMAN_HPP
# define MAXTREE_NAJMAN_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/sort_indexes.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/morpho/canonize.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      template <typename size_type>
      static
      inline
      size_type
      zfindroot(image2d<size_type>& parent, size_type p)
      {
        if (parent[p] == p)
          return p;
        else
          return parent[p] = zfindroot(parent, parent[p]);
      }

      template <typename size_type>
      struct aux_data
      {
	unsigned		rank;
	size_type	zpar;
      };

      template <typename size_type>
      static
      inline
      size_type
      zfindroot(image2d< aux_data<size_type> >& aux, size_type p)
      {
	size_type q = aux[p].zpar;
	if (p != q)
	  return aux[p].zpar = zfindroot(aux, q);
	else
	  return q;
      }

      template <typename size_type>
      static
      inline
      size_type
      mergeset(image2d< aux_data<size_type> >& aux, size_type p, size_type q)
      {
	if (aux[p].rank > aux[q].rank) std::swap(p,q);
	if (aux[p].rank == aux[q].rank) aux[q].rank += 1;
	aux[p].zpar = q;
	return q;
      }

    }

    template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
    std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
    maxtree_najman(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
    {
      typedef typename image2d<V>::size_type size_type;
      static constexpr size_type UNINITIALIZED = std::numeric_limits<size_type>::max();
      image2d<size_type>	parent;
      image2d<size_type>	lowestNode;
      image2d<internal::aux_data<size_type> >	        ttree;
      image2d<internal::aux_data<size_type> >	        tnode;

      resize(parent, ima).init(UNINITIALIZED);
      resize(lowestNode, ima);
      internal::aux_data<size_type> x = {0,0};
      resize(ttree, ima).init(x);
      resize(tnode, ima).init(x);

      /* coucou */

      std::vector<size_type> S = sort_indexes(ima, cmp);
      auto offsets = wrt_delta_index(ima, nbh.dpoints);
      std::cout << S.size() << std::endl;
      int j = S.size();
      for (int i = S.size()-1; i >= 0; --i)
	{
	  size_type p = S[i];
	  //std::cout << "Processing:" << p << " @ " << (int) ima[p] << std::endl;
	  // make set
	  {
	    parent[p] = p;
	    ttree[p].zpar = p;
	    tnode[p].zpar = p;
	    lowestNode[p] = p;
	  }


	  size_type curTree = internal::zfindroot(ttree, p); // zpar of p
	  size_type curNode = internal::zfindroot(tnode, lowestNode[curTree]); // zpar of p
	  for (unsigned k = 0; k < offsets.size(); ++k)
	    {
	      size_type q = p + offsets[k];
	      bool processed = (parent[q] != UNINITIALIZED);
	      if (processed)
		{
		  size_type adjTree = internal::zfindroot(ttree, q);
		  size_type adjNode = internal::zfindroot(tnode, lowestNode[adjTree]);
		  if (adjNode != curNode) // make union
		    {
		      //std::cout << "Mergin: " << root[r] << "->" << p << std::endl;
		      if (ima[curNode] == ima[adjNode]) {
			size_type tmp = internal::mergeset(tnode, adjNode, curNode);
			S[--j] = (tmp == curNode) ? adjNode : curNode;
			parent[adjNode] = parent[curNode] = tmp;
			curNode = tmp;
		      } else {
			parent[adjNode] = curNode;
			S[--j] = adjNode;
		      }
		      curTree = internal::mergeset(ttree, adjTree, curTree);
		      lowestNode[curTree] = curNode;
		    }
		}
	    }
	  //io::imprint(parent);
	}
      S[--j] = parent[S[0]];
      std::cout << j << std::endl;
      assert(j == 0);

      // canonization
      canonize(ima, S, parent);
      return std::make_pair(std::move(parent), std::move(S));
    }

  }

}

#endif // ! MAXTREE_NAJMAN_HPP
