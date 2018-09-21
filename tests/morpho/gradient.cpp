#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/win2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/morpho/structural/gradient.hpp>

#define BOOST_TEST_MODULE Morpho
#include <tests/test.hpp>

BOOST_AUTO_TEST_SUITE(gradient)

using namespace mln;

BOOST_AUTO_TEST_CASE(gradient_0)
{
  image2d<uint8> ima(10,10);
  iota(ima, 10);

  { // Fast: border wide enough
    auto win = make_rectangle2d(1, 3);
    auto out = morpho::structural::gradient(ima, win);

    static_assert( std::is_same<decltype(out)::value_type, int>::value,
                   "Error integral promotion should give int.");
    BOOST_CHECK(all(lor(out == 1, out == 2)));
  }
}

BOOST_AUTO_TEST_CASE(gradient_1)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  { // Fast: border wide enough
    auto win = make_rectangle2d(7, 7);
    auto out = morpho::structural::gradient(ima, win);
  }
}

// Border is not wide enough => call dilate + erode
BOOST_AUTO_TEST_CASE(gradient_2)
{
  image2d<uint8> ima(0);
  image2d<uint8> ima2;
  io::imread(MLN_IMG_PATH "small.pgm", ima);
  io::imread(MLN_IMG_PATH "small.pgm", ima2);

  auto win = make_rectangle2d(3, 3);
  auto out1 = morpho::structural::gradient(ima, win);
  auto out2 = morpho::structural::gradient(ima2, win);
  BOOST_CHECK(all(out1 == out2));
}


// Dilation on a with a vmorph / binary case
BOOST_AUTO_TEST_CASE(gradient_3)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);


  // Morpher has no extension
  auto win = make_rectangle2d(3, 3);
  auto out = morpho::structural::gradient(ima > 128, win);
}

// Dilation on a with a vmorph / binary case
BOOST_AUTO_TEST_CASE(gradient_4)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto win = make_rectangle2d(3, 3);
  image2d<uint8> out;
  resize(out, ima).init(0);
  auto tmp = out | where(ima > 128);
  morpho::structural::gradient(ima | where(ima > 128), win, std::less<uint8>(),
                   functional::l2norm_t<uint8>(), tmp);
}

 // On colors
BOOST_AUTO_TEST_CASE(gradient_5)
{
  image2d<rgb8> ima;
  image2d<rgb8> ima2(0);
  io::imread(MLN_IMG_PATH "small.ppm", ima);
  io::imread(MLN_IMG_PATH "small.ppm", ima2);

  auto win = make_rectangle2d(3, 3);
  auto out1 = morpho::structural::gradient(ima, win);
  auto out2 = morpho::structural::gradient(ima2, win);
  BOOST_CHECK(all(out1 == out2));
}

BOOST_AUTO_TEST_SUITE_END()

