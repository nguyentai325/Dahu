#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>

#include <mln/core/algorithm/equal.hpp>
#include <mln/core/algorithm/iota.hpp>


#define BOOST_TEST_MODULE Algorithm
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Equal)
{
using namespace mln;

image2d<uint8> ima(10, 10);
image2d<uint8> out(10, 10);

iota(ima, 0);
iota(out, 0);

BOOST_CHECK(equal(ima, out));

point2d p = {2,3};
out(p) = 12;
BOOST_CHECK(not equal(ima, out));
}
