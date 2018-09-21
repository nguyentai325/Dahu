#include <mln/core/image/image2d.hpp>
#include <mln/core/se/ball.hpp>
#include <mln/core/grays.hpp>
#include <mln/morpho/structural/dilate.hpp>
#include <mln/morpho/structural/erode.hpp>
#include <mln/io/imread.hpp>

#define BOOST_TEST_MODULE Morpho
#include <tests/test.hpp>

using namespace mln;

// Check that dilation/erosion are adjunctions
BOOST_AUTO_TEST_CASE(erode_0)
{
  image2d<uint8> ima;
  io::imread(MLN_IMG_PATH "small.pgm", ima);


  image2d<uint8> out1, out2;
  auto win = se::make_ball2d(3);

  auto comp = [](uint8 x) -> uint8 { return 255-x; };

  out1 = morpho::structural::dilate(ima, win);
  out2 = morpho::structural::erode(imtransform(ima, comp), win);

  MLN_CHECK_IMEQUAL(out2, imtransform(out1, comp));
}
