#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/equal.hpp>
#include <mln/morpho/extinction.hpp>


#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(extinction)

using namespace mln;

BOOST_AUTO_TEST_CASE(extinction_0)
{
  image2d<uint8> ima = { {0,3,2},
                         {2,2,3},
                         {2,3,1} };

  image2d<uint8> ref = { {3,0,1},
                         {0,0,0},
                         {0,0,2} };

  auto E = morpho::extinction(ima, c4);
  BOOST_CHECK(equal(E, ref));
}

BOOST_AUTO_TEST_CASE(extinction_1)
{
  image2d<uint8> ima = { {0,3,2},
                         {2,2,3},
                         {2,1,1} };

  image2d<uint8> ref = { {3,0,1},
                         {0,0,0},
                         {0,1,1} };

  auto E = morpho::extinction(ima, c4);
  BOOST_CHECK(equal(E, ref));
}



BOOST_AUTO_TEST_SUITE_END()

