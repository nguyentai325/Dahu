#ifndef APPS_SALIENCY_EXTINCTION_HPP
# define APPS_SALIENCY_EXTINCTION_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/trace.hpp>
# include <vector>

template <typename V, typename T, class Compare = std::less<V> >
mln::image2d<V>
extinction(const mln::image2d<V>& a,
	   const mln::image2d<T>& K,
	   const mln::image2d<unsigned>& parent,
	   const std::vector<unsigned>& S,
	   Compare cmp = Compare()
	   );


/************************************/
/**    Implementation		   **/
/************************************/

namespace internal
{

# ifndef INTERNAL_ZFINDROOT
#  define INTERNAL_ZFINDROOT

  unsigned
  zfindroot_(mln::image2d<unsigned>& par, unsigned x)
  {
    if (par[x] != x)
      par[x] = zfindroot_(par, par[x]);
    return par[x];
  }

# endif
}


template <typename V, typename T, class Compare>
mln::image2d<V>
extinction(const mln::image2d<V>& a,
	   const mln::image2d<T>& K,
	   const mln::image2d<unsigned>& parent,
	   const std::vector<unsigned>& S,
	   Compare cmp)
{
  using namespace mln;

  static const unsigned UNDEF = -1;

  trace::entering("extinction");

  // Retrieve the list of nodes
  // and sort them
  // also retrieve the child
  struct child_t {
    unsigned first_child = UNDEF;
    unsigned next_sibling = UNDEF;
  };

  std::vector<unsigned> nodes;
  image2d<child_t> childs;
  {
    resize(childs, K).init(child_t());

    nodes.reserve(S.size());
    nodes.push_back(S[0]);

    for (unsigned x: S)
      {
	unsigned q = parent[x];
	if (K[q] != K[x]) // Handle the root outside
	  {
	    // std::cout << "Setting node " << x << std::endl
	    // 	      << " parent: " << q << " / fc: " << childs[q].first_child << std::endl;
	    nodes.push_back(x);
	    if (childs[q].first_child == UNDEF)
	      childs[q].first_child = x;
	    else
	      {
		q = childs[q].first_child;
		while (childs[q].next_sibling != UNDEF)
		  q = childs[q].next_sibling;
		childs[q].next_sibling = x;
	      }
	  }
      }
    nodes.shrink_to_fit();

    std::sort(nodes.begin(), nodes.end(),
	      [&a, cmp](unsigned x, unsigned y) { return cmp(a[x], a[y]); });
  }


  // Compute extinction values.
  image2d<V> extmap;
  image2d<unsigned> par; // union-find struct

  resize(extmap, K).init(0);
  resize(par, K).init(UNDEF);


  for (unsigned x: nodes)
    {

      //std::cout << "Setting node" << x << std::endl;
      // Make set.
      unsigned zx = x;
      par[x] = zx;

      unsigned cchild = childs[x].first_child;
      // U-F
      // parent
      unsigned n = parent[x];
      //std::cout << "Doing: " << x << " > " << a[x] << std::endl;
      do
      {
	//std::cout << "  Nbh = " << n << " > " << a[n] << std::endl;
	if (par[n] != UNDEF) // deja_vu n ?
	  {
	    // zn is the local minimum of the neighboring component
	    unsigned zn = ::internal::zfindroot_(par, n);


	    // The highest minimum merges into the lowest minimum
	    // zn is the lowest
	    // zn in the highest (if not swap)
	    if (cmp(a[zx], a[zn]))
	      std::swap(zn, zx);

	    par[zx] = zn;
	    extmap[zx] = std::abs(a[x] - a[zx]);
	    //if (extmap[zx] != 0)
	    //  std::cout << "Setting minima " << zx << " : " << extmap[zx] << std::endl;
	    zx = zn;
	  }
      } while ( (n = cchild) != UNDEF and (cchild = childs[cchild].next_sibling, true) );
      // ^-- loop on children
    }
  // root extinction
  unsigned m = ::internal::zfindroot_(par, nodes[0]);
  extmap[m] = std::abs(a[nodes.back()] - a[m]);

  //std::cout << "Setting glob minimum " << m << " : " << extmap[m] << std::endl;

  trace::exiting();
  return extmap;
}




#endif // ! EXTINCTION_HPP
