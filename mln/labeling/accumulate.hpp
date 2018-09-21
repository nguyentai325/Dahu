#ifndef ACCUMULATE_HPP
# define ACCUMULATE_HPP

# include <mln/core/image/image.hpp>
# include <mln/accu/accumulator.hpp>
# include <type_traits>

namespace mln
{

  namespace labeling
  {

    template <typename I, typename Label, typename AccuLike>
    std::vector< typename accu::result_of<AccuLike, mln_point(I)>::type >
    p_accumulate(const Image<I>& lbl, Label nlabels, const AccumulatorLike<AccuLike>& accu_)
    {
      static_assert( std::is_same<mln_value(I), Label>::value,
		     "Image value type and Label type must match." );


      const I& ima = exact(lbl);
      auto acc = accu::make_accumulator(exact(accu_), mln_point(I) ());
      acc.init();

      typedef decltype(acc) Accu;
      std::vector<Accu> accumulators(nlabels+1, acc);

      // accumulate
      {
	mln_pixter(px, ima);
	mln_forall(px)
	  accumulators[px->val()].take(px->point());
      }


      // extract results
      typedef typename accu::result_of<AccuLike, mln_point(I)>::type R;
      std::vector<R> results(nlabels+1);
      {
	for (unsigned i = 0; i <= nlabels; ++i)
	  results[i] = accumulators[i].to_result();
      }

      return results;
    }

    template <typename I, typename J, typename Label, typename AccuLike>
    std::vector< typename accu::result_of<AccuLike, mln_value(J)>::type >
    v_accumulate(const Image<I>& lbl_,
		 const Image<J>& f_,
		 Label nlabels,
		 const AccumulatorLike<AccuLike>& accu_)
    {
      static_assert( std::is_same<mln_value(I), Label>::value,
		     "Image value type and Label type must match." );


      const I& lbl = exact(lbl_);
      const J& f = exact(f_);
      auto acc = accu::make_accumulator(exact(accu_), mln_value(J) ());
      acc.init();

      typedef decltype(acc) Accu;
      std::vector<Accu> accumulators(nlabels+1, acc);

      // accumulate
      {
	mln_pixter(px1, px2, lbl, f);
	mln_forall(px1, px2)
	  accumulators[px1->val()].take(px2->val());
      }

      // extract results
      typedef typename accu::result_of<AccuLike, mln_value(J)>::type R;
      std::vector<R> results(nlabels+1);
      {
	for (unsigned i = 0; i <= nlabels; ++i)
	  results[i] = accumulators[i].to_result();
      }

      return results;
    }



  }

}



#endif // ! ACCUMULATE_HPP
