#ifndef MLLN_CORE_CANVAS_ACCFPFN_HPP
# define MLLN_CORE_CANVAS_ACCFPFN_HPP

# include <mln/core/image/image.hpp>
# include <mln/accu/accumulator.hpp>
# include <type_traits>

namespace mln
{

  template <typename InputImage, typename OutputImage, typename Neighborhood, typename A, typename Extract>
  void
  accfpfn(const Image<InputImage>& ima,
	  const Neighborhood& nbh,
	  Image<OutputImage>&& out,
	  bool boundcheck = true,
	  Accumulator<A> acc = A(),
	  const Extract& extract = Extract());

  template <typename InputImage,
	    typename OutputImage,
	    typename Neighborhood,
	    typename A,
	    typename Extract>
  void
  accfpfn(const Image<InputImage>& ima,
	  const Neighborhood& nbh,
	  Image<OutputImage>& out,
	  bool boundcheck = true,
	  Accumulator<A> acc = A(),
	  const Extract& extract = Extract());

  template <typename InputImage,
	    typename Neighborhood,
	    typename A,
	    typename Extract>
  mln_ch_value(InputImage, typename std::result_of< Extract(A) >::type )
    accfpfn(const Image<InputImage>& ima,
	    const Neighborhood& nbh,
	    bool boundcheck = true,
	    Accumulator<A> acc = A(),
	    const Extract& extract = Extract());

  /************************/
  /** Implementation      */
  /************************/


  namespace impl
  {

    template <typename I, typename O, typename Neighborhood, typename A, typename Extract>
    void
    accfpfn(const I& ima_, const Neighborhood& nbh, O&& out_, bool boundcheck, A accu, const Extract& extract = Extract())
    {
      const I& ima = exact(ima_);
      O& out = exact(out_);

      mln_pixter(pin, pout, ima, out);
      mln_iter(qx, nbh(pin));

      if (boundcheck)
	{
	  mln_forall(pin, pout)
	    {
	      accu.init();
	      mln_forall(qx)
		if (ima.domain().has(qx->point()))
		  accu.take(qx->val());
	      pout->val() = extract(accu);
	    }
	}
      else
	{
	  mln_forall(pin, pout)
	    {
	      accu.init();
	      mln_forall(qx)
		accu.take(qx->val());
	      pout->val() = extract(accu);
	    }
	}
    }

  }

  template <typename InputImage, typename OutputImage, typename Neighborhood, typename A, typename Extract>
  void
  accfpfn(const Image<InputImage>& ima, const Neighborhood& nbh, Image<OutputImage>& out, bool boundcheck = true, Accumulator<A> acc = A(), const Extract& extract = Extract())
  {
    impl::accfpfn(exact(ima), nbh, exact(out), boundcheck, exact(acc), extract);
  }

  template <typename InputImage, typename OutputImage, typename Neighborhood, typename A, typename Extract>
  void
  accfpfn(const Image<InputImage>& ima, const Neighborhood& nbh, Image<OutputImage>&& out, bool boundcheck = true, Accumulator<A> acc = A(), const Extract& extract = Extract())
  {
    impl::accfpfn(exact(ima), nbh, exact(out), boundcheck, exact(acc), extract);
  }

  template <typename InputImage, typename Neighborhood, typename A, typename Extract>
  mln_ch_value(InputImage, typename std::result_of< Extract(A) >::type )
    accfpfn(const Image<InputImage>& ima, const Neighborhood& nbh, bool boundcheck = true, Accumulator<A> acc = A(), const Extract& extract = Extract())
  {
    mln_ch_value(InputImage, typename std::result_of< Extract(A) >::type )  out;
    resize(out, ima);
    accfpfn(ima, nbh, out, boundcheck, exact(acc), extract);
    return out;
  }

}


#endif // ! MLLN_CORE_CANVAS_ACCFPFN_HPP
