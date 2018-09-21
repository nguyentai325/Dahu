#ifndef OBJDETECTION_HPP
# define OBJDETECTION_HPP

# include <mln/core/image/image2d.hpp>
# include <vector>

namespace mln
{

  inline
  unsigned
  zfindroot(image2d<unsigned>& par, unsigned x)
  {
    if (par[x] != x)
      par[x] = zfindroot(par, par[x]);
    return par[x];
  }


  template <typename V>
  void
  getObjects(image2d<float>& energy,
	     image2d<V>& K,
	     image2d<unsigned>& parent,
	     const std::vector<unsigned>& S,
	     unsigned seuil)
  {
    static const unsigned NullNode = -1;

    struct child_t {
      unsigned first_child = -1;
      unsigned next_sibling = -1;
    };

    trace::entering("Meaningfull line detection");

    // Copy node indexes and set child relation
    std::vector<unsigned> nodes;
    image2d<child_t> crel;
    resize(crel, K);

    for (unsigned x: S)
      if (K[parent[x]] != K[x] or parent[x] == x) {
	nodes.push_back(x);
	unsigned y = crel[parent[x]].first_child;
	if (y == NullNode)
	  crel[parent[x]].first_child = x;
	else  {
	  while (crel[y].next_sibling != NullNode)
	    y = crel[y].next_sibling;
	  crel[y].next_sibling = x;
	}
      }

    // Sort
    std::sort(nodes.begin(), nodes.end(), [&energy](unsigned x, unsigned y) { return energy[x] < energy[y]; });



    // Union-find
    // par is union-find structure in *shape spaces*
    // the root of the "node component" is actually the
    // node with minimal energy
    static const unsigned UNSEEN = -1;
    image2d<unsigned> par, zpar;
    image2d<unsigned> area;
    image2d<unsigned> pmin;
    resize(zpar, K);
    resize(par, K).init(UNSEEN);
    resize(area, K);
    resize(pmin, K);

    unsigned num_node_deleted = 0;
    unsigned num_node_kept = 0;

    for (unsigned x: nodes)
      {
	par[x] = x;
	zpar[x] = x;
	area[x] = 1;
	pmin[x] = x;
	// Visist Nbh
	// Parent First
	{
	  unsigned q = parent[x];
	  if (par[q] != UNSEEN) {
	    unsigned y = zfindroot(zpar, q);
	    if (y != x) {
	      par[y] = x; // Caution: we choose y as the root (since energy[y] < energy[x])
	      zpar[y] = x;
	      area[x] += area[y];
	      pmin[x] = pmin[y];
	    }
	  }
	}
	// Children next
	{
	  unsigned q = crel[x].first_child;
	  if (q != NullNode)
	    do {
	      if (par[q] != UNSEEN) { //deja_vu
		unsigned y = zfindroot(zpar, q);
		if (y != x) {
		  if (energy[pmin[y]] < energy[pmin[x]])
		    pmin[x] = pmin[y];
		  par[y] = x;
		  zpar[y] = x;
		  area[x] += area[y];
		}
	      }
	      q = crel[q].next_sibling;
	    } while (q != NullNode);
	}
      }

    // Close energy
    {
      for (int i = nodes.size()-1; i >= 0; --i)
	{
	  unsigned x = nodes[i];
	  if (area[x] < seuil) {
	    energy[x] = energy[par[x]];
	    pmin[x] = pmin[par[x]];
	  }
	}
    }

    // simplify, remove non significant level lines
    // i.e keap only local minima whose
    // AND RECANONIZE !
    {
      for (unsigned x: S) {
	if (par[x] == UNSEEN) // non-canonical node
	  K[x] = K[parent[x]];
	else if (x != pmin[x]) { // not a local mimimum
	  K[x] = K[parent[x]];
	  energy[x] = energy[pmin[x]];
	  num_node_deleted++;
	} else {
	  num_node_kept++;
	  //std::cout << x << " ! " << energy[x] << std::endl;
	}

	unsigned q = parent[x];
	if (K[q] == K[parent[q]])
	  parent[x] = parent[q];
      }
    }
    std::cout << "==============" << std::endl;

    std::cout << "Level line removed: " << num_node_deleted << std::endl
	      << "Level line kept: " << num_node_kept << std::endl;

    trace::exiting();
  }

} // end of mln

#endif // ! OBJDETECTION_HPP
