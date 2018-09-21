#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>

#include <mln/core/algorithm/fill.hpp>


#define BOOST_TEST_MODULE Algorithm
#include <boost/test/unit_test.hpp>
#include <boost/range/algorithm/count.hpp>

BOOST_AUTO_TEST_CASE(Fill)
{
  using namespace mln;


  image2d<uint8> ima(10, 10);
  fill(ima, 69);

  mln_viter(v, ima);
  mln_forall(v)
    BOOST_CHECK(*v == 69);
}
