#include <array>
#include <mln/core/range/range.hpp>
#include <mln/core/forall.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>



BOOST_AUTO_TEST_CASE(stdrange_compatibility)
{

  using namespace mln;

  std::array<int, 4> x = {2, 5, 15, 22};

  BOOST_CHECK_EQUAL(rng::size(x), 4);

  {
    int i = 0;
    mln_foreach(int v, x)
      BOOST_CHECK_EQUAL(v, x[i++]);
  }

  {
    int i = 0;
    mln_iter(pv, x);
    mln_forall(pv)
      BOOST_CHECK_EQUAL(*pv, x[i++]);
  }

}
