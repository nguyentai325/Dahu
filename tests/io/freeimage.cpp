#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/grays.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/io/imprint.hpp>
#include <mln/core/algorithm/equal.hpp>
#include <mln/core/algorithm/iota.hpp>

#ifndef MLN_IMG_PATH
# error "MLN_IMG_PATH must be defined."
#endif


#define BOOST_TEST_MODULE IO
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE(Freeimage_pgm)
{
  using namespace mln;

  image2d<uint8> ima;
  image2d<uint8> ref(5, 5);

  iota(ref, 1);
  io::imread(MLN_IMG_PATH "/iota2d.pgm", ima);
  BOOST_CHECK(equal(ima, ref));
  io::imsave(ref, "test.tiff");
  io::imread("test.tiff", ima);
  BOOST_CHECK(equal(ima, ref));
}



BOOST_AUTO_TEST_CASE(Freeimage_ppm)
{
  using namespace mln;

  image2d<rgb8> ima;
  image2d<rgb8> ref(5, 5);

  mln_foreach(const image2d<rgb8>::pixel_type& pix, ref.pixels()) {
    pix.val()[0] = pix.point()[0];
    pix.val()[1] = pix.point()[1];
  }

  io::imread(MLN_IMG_PATH "/iota2d.ppm", ima);
  BOOST_CHECK(equal(ima, ref));
  io::imsave(ima, "test.tiff");
  io::imread("test.tiff", ima);
  BOOST_CHECK(equal(ima, ref));
}

BOOST_AUTO_TEST_CASE(Freeimage_pbm)
{
  using namespace mln;

  image2d<bool> ima;
  image2d<bool> ref(5, 5);

  mln_foreach(point2d p, ref.domain())
    ref(p) = ((p[0] % 2) == (p[1] % 2));

  io::imsave(ref, "test.tiff");
  io::imread("test.tiff", ima);
  BOOST_CHECK(equal(ima, ref));
}

BOOST_AUTO_TEST_CASE(Freeimage_slow_pgm)
{
  using namespace mln;

  image2d<uint8> ima;
  image2d<uint8> ref(5, 5);

  iota(ref, 1);
  io::imread(MLN_IMG_PATH "/iota2d.pgm", ima);
  BOOST_CHECK(equal(ima, ref));

  auto tmp = 2u * ref;
  io::imsave(imcast<uint32>(tmp), "test.tiff");
  image2d<unsigned> ima2;
  io::imread("test.tiff", ima2);
  BOOST_CHECK(equal(ima2, tmp));
}



BOOST_AUTO_TEST_CASE(Freeimage_slow_ppm)
{
  using namespace mln;

  image2d<rgb8> ima;
  image2d<rgb8> ref(5, 5);

  mln_foreach(const image2d<rgb8>::pixel_type& pix, ref.pixels()) {
    pix.val()[0] = pix.point()[0];
    pix.val()[1] = pix.point()[1];
  }

  io::imread(MLN_IMG_PATH "/iota2d.ppm", ima);
  BOOST_CHECK(equal(ima, ref));

  auto tmp = 2u * ref;
  io::imsave(imcast<rgb8>(tmp), "test.tiff");
  io::imread("test.tiff", ima);
  BOOST_CHECK(equal(ima, 2 * ref));
}

BOOST_AUTO_TEST_CASE(Freeimage_slow_pbm)
{
  using namespace mln;

  image2d<bool> ima;
  image2d<bool> ref(5, 5);

  mln_foreach(point2d p, ref.domain())
    ref(p) = ((p[0] % 2) == (p[1] % 2));

  io::imsave(lnot(ref), "test.tiff");
  io::imread("test.tiff", ima);
  BOOST_CHECK(equal(ima, lnot(ref)));
}

BOOST_AUTO_TEST_SUITE_END()
