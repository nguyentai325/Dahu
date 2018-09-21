#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/win2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/io/imread.hpp>
#include <mln/morpho/structural/opening.hpp>
#include <mln/morpho/structural/closing.hpp>

#define BOOST_TEST_MODULE Morpho
#include <tests/test.hpp>

BOOST_AUTO_TEST_SUITE(opening_closing)

using namespace mln;

BOOST_AUTO_TEST_CASE(opening_0)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto win = make_rectangle2d(3, 3);
  {
    auto out1 = morpho::structural::opening(ima, win);
    auto out2 = morpho::structural::closing(ima, win);

    BOOST_CHECK(all(out1 <= ima)); // anti-extensive
    BOOST_CHECK(all(out2 >= ima)); // extensive
  }
}

BOOST_AUTO_TEST_CASE(opening_1)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto comp = [](uint8 x) -> uint8 { return 255-x; };
  auto win = make_rectangle2d(3, 3);
  {
    auto out1 = morpho::structural::opening(imtransform(ima, comp), win);
    auto out2 = morpho::structural::closing(ima, win);

    MLN_CHECK_IMEQUAL(out1, imtransform(out2, comp)); // duality
  }
}

BOOST_AUTO_TEST_SUITE_END()

