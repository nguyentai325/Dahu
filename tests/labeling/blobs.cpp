#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/labeling/blobs.hpp>


#define BOOST_TEST_MODULE Labeling
#include <boost/test/unit_test.hpp>
#include <mln/io/imprint.hpp>

BOOST_AUTO_TEST_CASE(blobs_fast)
{
  using namespace mln;

  image2d<uint8> ima(5,5);
  iota(ima, 0);

  image2d<bool> mask = eval(ima % 2 == 0);

  image2d<uint8> lbl;
  unsigned nlabel;
  std::tie(lbl, nlabel) = labeling::blobs(mask, c4, uint8 ());
  BOOST_CHECK_EQUAL(nlabel, 13);
}


BOOST_AUTO_TEST_CASE(blobs_custom)
{

  using namespace mln;

  image2d<uint8> ima(5,5);
  iota(ima, 0);

  image2d<uint8> lbl;
  unsigned nlabel;

  std::tie(lbl, nlabel) = labeling::blobs(ima % 2 == 0, c4, uint8 ());
  BOOST_CHECK_EQUAL(nlabel, 13);

  std::tie(lbl, nlabel) = labeling::blobs(ima % 2 == 0, c8, uint8 ());
  BOOST_CHECK_EQUAL(nlabel, 1);
}
