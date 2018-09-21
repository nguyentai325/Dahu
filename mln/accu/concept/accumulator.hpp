#ifndef MLN_ACCU_CONCEPT_ACCUMULATOR_HPP
# define MLN_ACCU_CONCEPT_ACCUMULATOR_HPP

# include <mln/core/concept/object.hpp>
# include <mln/core/concept/check.hpp>
# include <boost/concept_check.hpp>

namespace mln
{

  /// \brief Concept for accumulator-like objects
  /// accumulator-like = accumulator | feature-set
  template <typename Acc>
  struct AccumulatorLike : Object<Acc>
  {
  protected:
    AccumulatorLike() = default;
    AccumulatorLike(const AccumulatorLike&) = default;
    AccumulatorLike& operator= (const AccumulatorLike&) = default;
  };

  /// \brief Concept for accumulator objects
  template <typename Acc>
  struct Accumulator : AccumulatorLike<Acc>
  {
    BOOST_CONCEPT_ASSERT((Accumulator<Acc>));
  };

  /// \brief Concept for feature-set objects
  template <typename E>
  struct FeatureSet : AccumulatorLike<E>
  {
  private:
    template <typename A>
    struct check_inner_struct_apply
    {
      template <typename T>
      using apply = typename A::template apply<T>::type;
    };

  public:
    BOOST_CONCEPT_USAGE(FeatureSet)
    {
      typedef typename E::features    features __attribute__((unused));
      check_inner_struct_apply<E> ();
    }
  };



  template <typename Acc>
  struct Accumulator_
  {
    typedef typename Acc::argument_type    argument_type;
    typedef typename Acc::result_type      result_type;


    BOOST_CONCEPT_USAGE(Accumulator_)
    {
      accu.init();
      accu.take(x);
      accu.take(accu2);
      res = accu.to_result();
    }

  private:
    Acc           accu, accu2;
    argument_type x;
    result_type   res;
  };

} // end of namespace mln

#endif //!MLN_ACCU_CONCEPT_ACCUMULATOR_HPP
