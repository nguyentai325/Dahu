#ifndef MAXTREE1D_HPP
# define MAXTREE1D_HPP

# include <mln/core/value/value_traits.hpp>
# include <mln/core/value/indexer.hpp>

namespace mln
{

  namespace morpho
  {

    template <typename V, typename StrictWeakOrdering>
    void
    maxtree1d(const image2d<V>& ima,
	      image2d<point2d>& parent, int row, StrictWeakOrdering cmp)
    {
      static const int nvalues = 1 << value_traits<V>::quant;
      point2d p = ima.domain().pmin;
      p[0] = row;

      int ncols = ima.ncols();
      point2d stack[nvalues];
      int sz = 0;
      point2d prec = p;
      ++p[1];
      for (int i = 1; i < ncols; ++i, ++p[1])
	{
	  //std::cout << "Processing: " << p << std::endl;
	  if (cmp(ima(prec),ima(p))) // ima(prec) < ima(p) => start new component
	    {
	      //std::cout << "  Push: " << prec << "in stack." << std::endl;
	      stack[sz++] = prec;
	      parent(prec) = prec;
	      prec = p;
	    }
	  else if (not cmp(ima(p), ima(prec))) // ima(p) == ima(prec) => extend component
	    {
	      //std::cout << "  Attach: " << p << " @ " << prec << std::endl;
	      parent(p) = prec;
	    }
	  else // ima(p) < ima(prec) => we need to attach prec to its parent
	    {
	      while (sz > 0 and not cmp(ima(stack[sz-1]), ima(p)))
		{
		  // std::cout << "  Pop: " <<  stack[sz-1] << " from stack." << std::endl;
		  // std::cout << "  Attach: " << prec << " @ " << stack[sz-1]  << std::endl;
		  parent(prec) = stack[sz-1];
		  prec = stack[sz-1];
		  --sz;
		}
	      // we have ima(p) <= ima(prec)
	      if (cmp(ima(p), ima(prec))) // ima(p) < ima(prec) => attach prec to p, p new component
		{
		  //std::cout << "  Attach: " << prec << " @ " << p << std::endl;
		  parent(prec) = p;
		  prec = p;
		}
	      else                        // ima(p) == ima(prec) => attach p to prec (canonization)
		{
		  //std::cout << "  Attach: " << p << " @ " << prec << std::endl;
		  parent(p) = prec;
		}
	    }
	}

      // Attach last point (i.e prec)
      while (sz > 0)
	{
	  // std::cout << "  Pop: " <<  stack[sz-1] << " from stack." << std::endl;
	  // std::cout << "  Attach: " << prec << " @ " << stack[sz-1]  << std::endl;
	  parent(prec) = stack[sz-1];
	  prec = stack[sz-1];
	  --sz;
	}
    }

  }

}



# ifndef MLN_INCLUDE_ONLY

# endif // ! MLN_INCLUDE_ONLY

#endif // ! MAXTREE1D_HPP
