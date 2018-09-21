#ifndef MLN_ACCU_ACCUMULATOR_BASE_HPP
# define MLN_ACCU_ACCUMULATOR_BASE_HPP

# include <mln/accu/accumulator.hpp>

namespace mln
{

  namespace accu
  {


    // \brief A base class for adaptable accumulator
    template <typename E, typename ArgumentType, typename ResultType, typename Feature>
    struct accumulator_base;





    /*********************/
    /** Implementation   */
    /********************/

    template <typename E, typename ArgumentType, typename ResultType, typename Feature>
    struct accumulator_base : Accumulator<E>
    {
      typedef ArgumentType argument_type;
      typedef ResultType   result_type;
      typedef Feature	   feature;
      typedef boost::mpl::set<Feature> provides;

      typedef std::false_type          has_untake;


      // You mut implement the following function
      // void init()
      // void take(const argument_type&)
      // void take(const Accumulator<A>& other)
      // friend result_type extract(const E&, Feature)
      //
      // and the optional following methods:
      // bool stop()
      // void untake(const argument_type&)

      constexpr bool stop() const
      {
        return false;
      }


      result_type to_result() const
      {
        return extract(exact(*this), feature ());
      }

    };

  }

}

#endif // ! MLN_ACCU_ACCUMULATOR_BASE_HPP
