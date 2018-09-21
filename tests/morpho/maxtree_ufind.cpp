#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>

#include <mln/core/algorithm/iota.hpp>
#include <mln/morpho/maxtree_ufind.hpp>
#include <mln/io/imprint.hpp>

#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>

#include <mln/core/value/int.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(maxtree_ufind)
{
  using namespace mln;

  image2d<uint8> ima(5,5);
  iota(ima, 2);

  io::imprint(ima);
  {
    image2d<std::size_t> parent = morpho::maxtree(ima, c4, std::less<uint8> ());
    io::imprint(parent);
  }

  {
    image2d<std::size_t> parent = morpho::maxtree(ima, c4, std::greater<uint8> ());
    io::imprint(parent);
  }

}
