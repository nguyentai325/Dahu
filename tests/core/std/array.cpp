
#include <mln/core/std/array.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Array)
{
  using namespace mln;

  array<int, 3> x = {2,3,4};

  BOOST_CHECK_EQUAL(x[0], 2);
}
