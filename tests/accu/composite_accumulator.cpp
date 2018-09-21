#include <mln/accu/composite_accumulator.hpp>
#include <mln/accu/accumulator.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>
#include <boost/mpl/has_key.hpp>

namespace mln
{

  namespace accu
  {
    namespace accumulators
    {
      template <typename T> struct accu1;
      template <typename T> struct accu2;
      template <typename T> struct accu3;
      template <typename T> struct accu4;
      template <typename T> struct accu5;
    }

    namespace features
    {

      struct f1 : simple_feature_facade<f1,  accumulators::accu1> {};
      struct f2 : simple_feature_facade<f2,  accumulators::accu2> {};
      struct f3 : simple_feature_facade<f3,  accumulators::accu3> {};
      struct f4 : simple_feature_facade<f4,  accumulators::accu4> {};
      struct f5 : simple_feature_facade<f5,  accumulators::accu5> {};


      template<>
      struct depends<f1> { typedef boost::mpl::set<f2> type; };

      template<>
      struct depends<f2> { typedef boost::mpl::set<f2,f3> type; };

      template<>
      struct depends<f3> { typedef boost::mpl::set<f4,f5> type; };

      template<>
      struct depends<f5> { typedef boost::mpl::set<f4> type; };

    }

    namespace accumulators
    {

#define DEF_ACCU(N)							\
      template <typename T>						\
      struct accu##N : composite_accumulator_facade< accu##N<T>, T, T, features::f##N> \
      {									\
	typedef T argument_type;					\
	typedef boost::mpl::set<features::f##N> provides;		\
	void init() {};							\
	void take(const T&) {};						\
	template <typename Other> void take(const Accumulator<Other>&) {} \
	friend T extract(const accu##N&, features::f##N) {return T(); };	\
      };


      DEF_ACCU(1);
      DEF_ACCU(2);
      DEF_ACCU(3);
      DEF_ACCU(4);
      DEF_ACCU(5);
    }

  }

}

BOOST_AUTO_TEST_CASE(Accu_Features)
{
  using namespace mln::accu;

  typedef boost::mpl::set<features::f1> Set;
  typedef typename internal::resolve_dependances<Set>::type Features;



  static_assert(boost::mpl::has_key<Features, features::f1>::value, "");
  static_assert(boost::mpl::has_key<Features, features::f2>::value, "");
  static_assert(boost::mpl::has_key<Features, features::f3>::value, "");
  static_assert(boost::mpl::has_key<Features, features::f4>::value, "");
  static_assert(boost::mpl::has_key<Features, features::f5>::value, "");


  auto x = make_accumulator(features::f3 (), int() );
  extract(x, features::f3 ());
  extract(x, features::f4 ());
  extract(x, features::f5 ());
}
