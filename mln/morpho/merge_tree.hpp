#ifndef MLN_CORE_MORPHO_MERGE_TREE_HPP
# define MLN_CORE_MORPHO_MERGE_TREE_HPP

# include <tbb/compat/thread>
# include <mln/morpho/maxtree_routines.hpp>

namespace mln
{

  namespace morpho
  {


    // tree merging
    template <typename V, typename StrictWeakOrdering>
    void merge_tree(const image2d<V>& ima,
		    image2d<typename image2d<V>::size_type>& parent,
		    box2d domain,
		    StrictWeakOrdering cmp)
    {
      typedef typename image2d<V>::size_type size_type;
      mln_precondition(!domain.empty());

      point2d p_ = domain.pmin;
      point2d q_ = domain.pmin;
      p_[0] = domain.pmax[0]-1;
      q_[0] = domain.pmax[0];
      size_type p = ima.index_of_point(p_);
      size_type q = ima.index_of_point(q_);

      unsigned ncols = ima.ncols();
      for (unsigned i = 0; i < ncols; ++i, ++p, ++q)
	{
	  size_type x = internal::zfind_repr(ima, parent, p);
	  size_type y = internal::zfind_repr(ima, parent, q);
	  if (cmp(ima[x], ima[y]))
	    std::swap(x, y);

	  while (x != y)
	    {
	      //std::cout << "-- Merge: " << x << " @ " << y << std::endl;
	      // check that x and y are representative
	      mln_assertion(x == parent[x] or ima[parent[x]] != ima[x]);
	      mln_assertion(y == parent[y] or ima[parent[y]] != ima[y]);
	      mln_assertion(!cmp(ima[x], ima[y])); // ima(y) <= ima(x)

	      // we want to attach x to y
	      if (parent[x] == x)
		{
		  parent[x] = y;
		  break;
		}
	      else
		{
		  size_type z = internal::zfind_parent(ima, parent, x);
		  if (!cmp(ima[z], ima[y])) // ima(y) <= ima(z)
		    x = z;
		  else
		    {
		      parent[x] = y;
		      x = y;
		      y = z;
		    }
		}
	    }

	}

    }

    template <typename size_type>
    bool
    check_S(const image2d<size_type>& parent, const size_type* begin, const size_type* end)
    {
      image2d<bool> dejavu;
      resize(dejavu, parent).init(false);

      dejavu[*begin] = true;
      for (;begin != end; ++begin) {
        assert(dejavu[parent[*begin]]);
	if (!dejavu[parent[*begin]])
	  return false;
        dejavu[*begin] = true;
      }
      return true;
    }

    /*
    static
    inline
    void
    addP(image2d<bool>& dejavu, const image2d<std::size_t>& parent, std::size_t p, (std::size_t*)& dst)
    {
      if (!deja_vu[p])
        {
          deja_vu[p] = true;
          addAncestor(deja_vu, parent, parent[p], dst);
          *(dst++) = p;
        }
    }


    // tree merging
    // note: dst may alias with src2 but not src1
    void merge_S(const image2d<std::size_t>& parent,
		 image2d<bool>& deja_vu,
		 const std::size_t* src1, const std::size_t* end1,
		 const std::size_t* src2, const std::size_t* end2,
		 std::size_t* dst)
    {
      mln_precondition(dst != src1);
      mln_precondition(src2 == end1 or (src2 == (dst + (end1-src1))));

      std::size_t* buffer = dst;
      std::size_t n = (end1-src1) + (end2-src2);

      if (parent[*src1] == *src1)
	*dst = *(src1++);
      else
	*dst = *(src2++);

      mln_assertion(parent[*dst] == *dst);
      deja_vu[*(dst++)] = true;

      for (;src1 != end1; ++src1)
        if (!deja_vu[*src1])
          {
            std::size_t q = parent[*src1];
            while (!deja_vu[q]) {
              *ds
            }
            *(dst++) = *src1;
          }

	  if (deja_vu[parent[*src1]])
	    *dst = *(src1++);
	  else if (deja_vu[parent[*src2]])
	    *dst = *(src2++);
          else
            {
              
            }
	  deja_vu[*(dst++)] = true;
	}

      mln_assertion((src1 == end1) != (src2 == end2));
      if (src1 != end1)
	std::copy(src1, end1, dst);
      else if (src2 != end2 and src2 != dst)
	std::copy(src2, end2, dst);

      check_S(parent, buffer, buffer+n);
    }
    */
  }

}


#endif // !MLN_CORE_MORPHO_MERGE_TREE_HPP
