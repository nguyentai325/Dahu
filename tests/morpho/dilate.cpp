#include <mln/core/image/image2d.hpp>
#include <mln/core/win2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/morpho/structural/dilate.hpp>
#include <mln/io/imread.hpp>



#define BOOST_TEST_MODULE Morpho
#include <tests/test.hpp>

BOOST_AUTO_TEST_SUITE(dilation)

using namespace mln;

BOOST_AUTO_TEST_CASE(dilate_0)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  { // Fast: border wide enough
    auto win = make_rectangle2d(7, 7);
    auto out = morpho::structural::dilate(ima, win);
    BOOST_CHECK( all(out >= ima) ); // extensive
  }
}

// Border is not wide enough => use morpher for bound checking
BOOST_AUTO_TEST_CASE(dilate_1)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto win = make_rectangle2d(9, 9);
  auto out = morpho::structural::dilate(ima, win);
  BOOST_CHECK( all(out >= ima) ); // extensive
}


// Dilation on a with a vmorph / binary case
BOOST_AUTO_TEST_CASE(dilate_2)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);


  auto win = make_rectangle2d(3, 3);
  auto out = morpho::structural::dilate(ima > 128, win);
  BOOST_CHECK( all(out >= (ima > 128) )); // extensive
}

// Dilation on a with a vmorph / binary case
BOOST_AUTO_TEST_CASE(dilate_view)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto win = make_rectangle2d(3, 3);
  auto out = clone(ima);
  auto tmp = out | where(ima > 128);
  morpho::structural::dilate(ima | where(ima > 128), win, tmp, std::less<uint8>() );
  BOOST_CHECK( all((out | (ima <= 128)) == (ima | (ima <= 128))) );
  BOOST_CHECK( all(out >= ima) ); // extensive
}

 // Custom comparison function, erosion
BOOST_AUTO_TEST_CASE(dilate_with_custom_cmp_function)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);

  auto win = make_rectangle2d(5, 5);
  auto out = morpho::structural::dilate(ima, win, std::greater<uint8> ());
  BOOST_CHECK( all(out <= ima) ); // anti-extensive
}


// Dilation of a binary image
BOOST_AUTO_TEST_CASE(dilate_binary)
{
  image2d<bool> ima(11,11);
  fill(ima, false);
  ima.at(0,0) = true;
  ima.at(5,5) = true;
  ima.at(10,10) = true;

  auto win = make_rectangle2d(3, 3);
  auto out = morpho::structural::dilate(ima, win);
  BOOST_CHECK( all(ima <= out) ); // anti-extensive
}

// Dilation of a bianry image
BOOST_AUTO_TEST_CASE(dilate_binary_2)
{
  image2d<bool> ima;
  io::imread(MLN_IMG_PATH "tiny.pbm", ima);

  auto win = make_rectangle2d(3, 3);
  auto out = morpho::structural::dilate(ima, win);
  BOOST_CHECK( all(ima <= out) ); // anti-extensive
}


// Dilation of a rgb image
BOOST_AUTO_TEST_CASE(dilate_rgb)
{
  image2d<rgb8> ima;
  io::imread(MLN_IMG_PATH "small.ppm", ima);

  auto win = make_rectangle2d(5, 5);
  auto out = morpho::structural::dilate(ima, win);
  BOOST_CHECK( all(red(ima) <= red(out)) ); // anti-extensive
  BOOST_CHECK( all(green(ima) <= green(out)) ); // anti-extensive
  BOOST_CHECK( all(blue(ima) <= blue(out)) ); // anti-extensive
}


 BOOST_AUTO_TEST_SUITE_END()

