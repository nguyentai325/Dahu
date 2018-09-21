#ifndef MLN_CORE_ALGORITHM_ACCUMULATE_HPP
# define MLN_CORE_ALGORITHM_ACCUMULATE_HPP

# include <mln/core/image/image.hpp>
# include <mln/accu/accumulator.hpp>

namespace mln {

  template <typename I, class AccuLike, class Extractor = accu::default_extractor>
  typename accu::result_of<AccuLike, mln_value(I), Extractor>::type
  accumulate(const Image<I>& input, const AccumulatorLike<AccuLike>& accu, const Extractor& ex = Extractor ());


  template <typename I, class BinaryOperator, class V>
  typename std::enable_if< !is_a<BinaryOperator, AccumulatorLike>::value, V>::type
  accumulate(const Image<I>& input, const BinaryOperator& op, V init);


  /*********************/
  /*** Implementation  */
  /*********************/

  template <typename I, class AccuLike, class Extractor>
  typename accu::result_of<AccuLike, mln_value(I), Extractor>::type
  accumulate(const Image<I>& input, const AccumulatorLike<AccuLike>& accu_, const Extractor& ex)
  {
    const I& ima = exact(input);
    auto a = accu::make_accumulator(exact(accu_), mln_value(I) ());

    mln_foreach(const auto& v, ima.values())
      a.take(v);

    return ex(a);
  }

  template <typename I, class BinaryOperator, class V>
  typename std::enable_if< !is_a<BinaryOperator, AccumulatorLike>::value, V>::type
  accumulate(const Image<I>& input, const BinaryOperator& op, V init)
  {
    const I& ima = exact(input);

    mln_foreach(const auto& v, ima.values())
      init = op(init, v);

    return init;
  }

}


#endif // ! MLN_CORE_ALGORITHM_ACCUMULATE_HPP
