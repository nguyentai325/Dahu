#include <mln/core/vec.hpp>
#include <mln/core/colors.hpp>
#include <mln/accu/accumulators/infsup.hpp>
#include <mln/accu/accumulators/minmax.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Infsup)
{
  using namespace mln;
  using namespace mln::accu;

  typedef vec3i Vec;
  typedef productorder_less<Vec> Compare;

  // product order comparison
  {
    accumulators::infsup<Vec> acc;

    acc.take(Vec(2,-5,6));
    acc.take(Vec(4,-1,-3));
    acc.take(Vec(-2,-5,7));

    BOOST_CHECK_EQUAL( extractor::inf(acc), Vec(-2,-5,-3));
    BOOST_CHECK_EQUAL( extractor::sup(acc), Vec(4,-1,7));
  }

  // defaut comparison on generic vectors is lexicographical
  {
    accumulators::minmax<Vec> acc;

    acc.take(Vec(4,-5,6));
    acc.take(Vec(4,-1,-3));
    acc.take(Vec(-2,-5,7));

    BOOST_CHECK_EQUAL( extractor::min(acc), Vec(-2,-5,7));
    BOOST_CHECK_EQUAL( extractor::max(acc), Vec(4,-1,-3));
  }

  // defaut comparison on generic vectors is lexicographical
  // This must not compile and generate an error message
  /*
  {
    typedef rgb<int> Vec;
    accumulators::infsup<Vec> acc;

    acc.take(Vec(4,-5,6));
    acc.take(Vec(4,-1,-3));
    acc.take(Vec(-2,-5,7));

    BOOST_CHECK_EQUAL( extractor::inf(acc), Vec(-2,-5,7));
    BOOST_CHECK_EQUAL( extractor::sup(acc), Vec(4,-1,-3));
  }
  */
}
