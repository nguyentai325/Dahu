#include <mln/core/vec.hpp>
#include <mln/core/vec/vec_math_ops.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>



BOOST_AUTO_TEST_CASE(vec_math_ops)
{
  using namespace mln;

  {
    vec3b v = {-127,2,-1};

    BOOST_CHECK_EQUAL(l0norm(v),   1);
    BOOST_CHECK_EQUAL(linfnorm(v), 127);
    BOOST_CHECK_EQUAL(l1norm(v), 130); // check type promotion
    BOOST_CHECK_EQUAL(l2norm(v), std::sqrt(127 * 127 + 4 + 1));
  }
}
