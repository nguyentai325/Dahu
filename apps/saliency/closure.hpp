#ifndef MLN_APPS_SALIENCY_CLOSURE_HPP
# define MLN_APPS_SALIENCY_CLOSURE_HPP

#include <mln/core/image/image2d.hpp>



template <typename V, class Compare = std::less<unsigned> >
mln::image2d<float>
area_close(const mln::image2d<float>& attr,
	   const mln::image2d<V>& K,
	   const mln::image2d<unsigned>& parent,
	   const std::vector<unsigned>& S,
	   unsigned lambda,
	   Compare cmp = Compare ());


template <typename V, class Compare = std::less<float> >
mln::image2d<float>
height_close(const mln::image2d<float>& attr,
	     const mln::image2d<V>& K,
	     const mln::image2d<unsigned>& parent,
	     const std::vector<unsigned>& S,
	     float lambda,
	     Compare cmp = Compare ());

/************************/
/** Implementation     **/
/************************/

namespace internal
{

  struct child_t
  {
    static constexpr unsigned UNDEF = -1;

    unsigned first_child = UNDEF;
    unsigned next_sibling = UNDEF;
  };


  template <typename V>
  mln::image2d<child_t>
  getchilds(const mln::image2d<V>& K,
	    const mln::image2d<unsigned>& parent,
	    const std::vector<unsigned>& S)
  {
    using namespace mln;
    image2d<child_t> childs;
    resize(childs, K).init(child_t());

    for (unsigned x: S)
      {
	unsigned q = parent[x];
	if (K[q] != K[x]) // Handle the root outside
	  {
	    if (childs[q].first_child == child_t::UNDEF)
	      childs[q].first_child = x;
	    else
	      {
		q = childs[q].first_child;
		while (childs[q].next_sibling != child_t::UNDEF)
		  q = childs[q].next_sibling;
		childs[q].next_sibling = x;
	      }
	  }
      }
    return childs;
  }

  template <typename V>
  std::vector<unsigned>
  getnodes(const mln::image2d<V>& K,
	   const mln::image2d<unsigned>& parent,
	   const std::vector<unsigned>& S)
  {
    std::vector<unsigned> nodes;
    nodes.reserve(S.size());

    nodes.push_back(S[0]);
    for (unsigned x: S)
      {
	unsigned y = parent[x];
	if (K[x] != K[y]) // x is a node
	  nodes.push_back(x);
      }
    nodes.shrink_to_fit();

    return nodes;
  }

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


template <typename V, class Compare>
mln::image2d<float>
area_close(const mln::image2d<float>& attr,
	   const mln::image2d<V>& K,
	   const mln::image2d<unsigned>& parent,
	   const std::vector<unsigned>& S,
	   unsigned lambda,
	   Compare cmp)
{
  using namespace mln;

  image2d< ::internal::child_t > childs = ::internal::getchilds(K, parent, S);
  std::vector<unsigned> nodes = ::internal::getnodes(K, parent, S);

  std::sort(nodes.begin(), nodes.end(), [cmp,&attr](unsigned x, unsigned y) {
      return cmp(attr[x], attr[y]);
    });

  // perform the closure.
  static const unsigned UNDEF = ::internal::child_t::UNDEF;

  image2d<unsigned> par, area;
  resize(par, K).init(UNDEF);
  resize(area, K);

  for (unsigned x: nodes)
    {
      // make set
      par[x] = x;
      area[x] = 1;

      unsigned cchild = childs[x].first_child;
      unsigned n = parent[x];

      do // foreach neighbor n in tree topo
	{
	  if (par[n] != UNDEF) // dejavu(n) ?
	    {
	      unsigned r = ::internal::zfindroot_(par, n);
	      if (r != x) // merge set
		{
		  if (area[r] < lambda) // r should be filtered
		    {
		      par[r] = x;
		      area[x] += area[r];
		    }
		  else
		    {
		      // r should merge with x but no need since
		      // r is not filtered. Just set Area[x] = lambda
		      // to say that x should not be filtered
		      area[x] = lambda;
		    }
		}
	    }
	}
      while ((n = cchild) != UNDEF and (cchild = childs[cchild].next_sibling, true) );
    }

  // filter now:
  image2d<float> out;
  resize(out, attr);
  out[nodes[nodes.size()-1]] = value_traits<unsigned, Compare>::max();
  for (int i = nodes.size()-1; i >= 0; --i)
    {
      unsigned x = nodes[i];
      if (area[x] < lambda)
	out[x] = out[par[x]];
      else
	out[x] = attr[x];
    }

  return out;
}


template <typename V, class Compare>
mln::image2d<float>
height_close(const mln::image2d<float>& attr,
	     const mln::image2d<V>& K,
	     const mln::image2d<unsigned>& parent,
	     const std::vector<unsigned>& S,
	     float lambda,
	     Compare cmp)
{
  using namespace mln;

  image2d< ::internal::child_t > childs = ::internal::getchilds(K, parent, S);
  std::vector<unsigned> nodes = ::internal::getnodes(K, parent, S);

  std::sort(nodes.begin(), nodes.end(), [cmp,&attr](unsigned x, unsigned y) {
      return cmp(attr[x], attr[y]);
    });

  // perform the closure.
  static const unsigned UNDEF = ::internal::child_t::UNDEF;

  image2d<unsigned> par;
  image2d<float> minimum;
  resize(par, K).init(UNDEF);
  resize(minimum, K);

  for (unsigned x: nodes)
    {
      // make set
      par[x] = x;
      minimum[x] = attr[x];

      unsigned cchild = childs[x].first_child;
      unsigned n = parent[x];

      float clevel = attr[x];

      do // foreach neighbor n in tree topo
	{
	  if (par[n] != UNDEF) // dejavu(n) ?
	    {
	      unsigned r = ::internal::zfindroot_(par, n);
	      if (r != x) // merge set
		{
		  if (std::abs(clevel - minimum[r]) < lambda) // r should be filtered
		    {
		      par[r] = x;
		      minimum[x] = std::min(minimum[x], minimum[r], cmp);
		    }
		  else
		    {
		      // r should merge with x but no need since
		      // r is not filtered. Just set Min[x] = ^(float, <)
		      // to say that x should not be filtered
		      minimum[x] = value_traits<float, Compare>::min();
		    }
		}
	    }
	}
      while ((n = cchild) != UNDEF and (cchild = childs[cchild].next_sibling, true) );
    }

  // filter now:
  image2d<float> out;
  resize(out, attr);
  out[nodes[nodes.size()-1]] = value_traits<float, Compare>::max();
  for (int i = nodes.size()-1; i >= 0; --i)
    {
      unsigned x = nodes[i];
      if ( std::abs(attr[x] - minimum[x]) < lambda )
	out[x] = out[par[x]];
      else
	out[x] = attr[x];
    }

  return out;
}





#endif // ! MLN_APPS_SALIENCY_CLOSURE_HPP
