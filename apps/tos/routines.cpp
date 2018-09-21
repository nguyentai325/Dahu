#include "routines.hpp"

unsigned
lca(const mln::image2d<unsigned>& depth,
    const mln::image2d<unsigned>& parent,
    unsigned x, unsigned y)
{
  if (x == -1)
    return y;
  else if (y == -1)
    return x;
  else
    {
      while (x != y) {
	if (depth[x] > depth[y])
	  x = parent[x];
	else if (depth[y] > depth[x])
	  y = parent[y];
	else {
	  x = parent[x];
	  y = parent[y];
	}
      }
      return x;
    }
}
