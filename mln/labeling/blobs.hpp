#ifndef MLN_LABELING_BLOBS_HPP
# define MLN_LABELING_BLOBS_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/trace.hpp>
# include <vector>
# include <type_traits>

namespace mln
{

  namespace labeling
  {

    /// \brief labelize connected components of a binary image ima.
    template <typename I, typename Neighborhood, typename Label = unsigned>
    std::pair< mln_ch_value(I, Label), Label >
    blobs(const Image<I>& ima, Neighborhood nbh, Label lbl = Label());


    /******************************/
    /*** Implementation         ***/
    /******************************/

    namespace impl
    {

      namespace generic
      {

	template <typename I, typename Neighborhood, typename Label, typename O>
	Label
	blobs_no_boundcheck(const I& ima, Neighborhood nbh, Label lbl, O& out)
	{
	  typedef mln_value(I) V;
	  typedef mln_point(I) P;

	  mln_entering("mln::labeling::impl::generic::blobs_no_boundcheck");

	  Label bg = lbl;
	  std::vector<P> queue;
	  queue.reserve(ima.domain().size());

	  extension::fill(out, value_traits<Label>::max());

	  P q;
	  mln_iter(n, nbh(q));

	  mln_foreach(P p, ima.domain())
	    {
	      if (ima(p) and out(p) == bg)
		{
		  queue.push_back(p);
		  ++lbl;
		  mln_assertion(lbl <= value_traits<Label>::max());
		  while (not queue.empty())
		    {
		      q = queue.back();
		      queue.pop_back();
		      out(q) = lbl;
		      mln_forall(n)
			if (out.at(*n) == bg and ima(*n))
			  queue.push_back(*n);
		    }
		}
	    }

	  mln_exiting();

	  return lbl;
	}

	template <typename I, typename Neighborhood, typename Label, typename O>
	Label
	blobs_boundcheck(const I& ima, Neighborhood nbh, Label lbl, O& out)
	{
	  typedef mln_value(I) V;
	  typedef mln_point(I) P;

	  mln_entering("mln::labeling::impl::generic::blobs_boundcheck");

	  Label bg = lbl;

	  std::vector<P> queue;
	  queue.reserve(ima.domain().size());

	  P q;
	  mln_iter(n, nbh(q));

	  mln_foreach(P p, ima.domain())
	    {
	      if (ima(p) and out(p) == bg)
		{
		  queue.push_back(p);
		  ++lbl;
		  mln_assertion(lbl <= value_traits<Label>::max());
		  while (not queue.empty())
		    {
		      q = queue.back();
		      queue.pop_back();
		      out(q) = lbl;
		      mln_forall(n)
			if (ima.domain().has(*n) and ima(*n) and out(*n) == bg)
			  queue.push_back(*n);
		    }
		}
	    }

	  mln_exiting();

	  return lbl;
	}

      }

    }


    template <typename I, typename Neighborhood, typename Label>
    std::pair< mln_ch_value(I, Label), Label >
    blobs(const Image<I>& ima_, Neighborhood nbh, Label lbl)
    {
      typedef mln_value(I) V;
      typedef mln_point(I) P;
      static_assert(std::is_same<V, bool>::value, "Only supports binary image (type: bool)");

      mln_entering("mln::labeling::blobs");

      const I& ima = exact(ima_);
      Label bg = lbl;

      int status;
      mln_ch_value(I, Label) out = imchvalue<Label>(ima).adjust(nbh).init(bg).get_status(status);

      if (status == 0)
        lbl = impl::generic::blobs_no_boundcheck(ima, nbh, lbl, out);
      else
        lbl = impl::generic::blobs_boundcheck(ima, nbh, lbl, out);

      mln_exiting();

      return std::make_pair(out, lbl);
    }

  }

}

#endif // ! MLN_LABELING_BLOBS_HPP
