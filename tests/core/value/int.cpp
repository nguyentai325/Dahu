#include <mln/core/value/int.hpp>
#include <type_traits>

#define BOOST_TEST_MODULE Value
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Integers)
{
  using namespace mln;

  UInt<12> x = 13;
  UInt<12> y = 18;
  BOOST_CHECK_EQUAL(x + y, 31);

  Int<14> z = 12;
  y += z;
  ++y;
  BOOST_CHECK_EQUAL(y, 31);

  Int<7> zz = x + z;
  BOOST_CHECK_EQUAL(zz, 25);

  typedef typename std::common_type<UInt<12>, UInt<12> >::type T;
  BOOST_CHECK((std::is_same<T, UInt<12>>::value));
  BOOST_CHECK((sizeof(UInt<12>) == sizeof(short)));
}


