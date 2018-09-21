#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>

#include <mln/core/algorithm/clone.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/extension/fill.hpp>

#include <mln/io/imprint.hpp>

#define BOOST_TEST_MODULE Extension
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Fill)
{
  using namespace mln;


  image2d<uint8> ima(10, 10);
  iota(ima, 0);
  image2d<uint8> out = clone(ima);

  extension::fill(ima, 42);

  io::imprint_with_border(ima);
  BOOST_CHECK( all(ima == out) );
}
