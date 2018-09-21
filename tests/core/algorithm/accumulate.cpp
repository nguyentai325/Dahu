#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/accu/accumulators/sum.hpp>
#include <mln/accu/accumulators/min.hpp>
#include <mln/accu/accumulators/max.hpp>

#define BOOST_TEST_MODULE Algorithm
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(Accumulate_1)
{
  using namespace mln;

  image2d<uint8> ima(10, 10);
  iota(ima, 0);

  // Expected overflow
  {
    int res = accumulate(ima, std::plus<uint8>(), 0);
    BOOST_CHECK_EQUAL(res, ((99*100) / 2) % 256);
  }

  // No overflow
  {
    int res = accumulate(ima, std::plus<int>(), 0);
    BOOST_CHECK_EQUAL(res, ((99*100) / 2));
  }
}


BOOST_AUTO_TEST_CASE(Accumulate_2)
{
  using namespace mln;

  image2d<uint8> ima(10, 10);
  iota(ima, 0);

  // No overflow (uint8 + uint8 -> int)
  {
    int res = accumulate(ima, accu::features::sum<> ());
    BOOST_CHECK_EQUAL(res, ((99*100) / 2));
  }
}


BOOST_AUTO_TEST_CASE(Accumulate_3)
{
  using namespace mln;

  image2d<uint8> ima(10, 10);
  iota(ima, 0);

  // No overflow (uint8 + uint8 -> int)
  {
    auto acc = accumulate(ima, accu::features::min<> () & accu::features::max<> ());
    BOOST_CHECK_EQUAL(accu::extractor::min(acc), 0);
    BOOST_CHECK_EQUAL(accu::extractor::max(acc), 99);
  }
}
