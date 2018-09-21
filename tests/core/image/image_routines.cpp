#include <mln/core/image/image2d.hpp>
#include <mln/core/algorithm/iota.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>
#include <mln/io/imprint.hpp>

BOOST_AUTO_TEST_CASE(Where_binary)
{
  using namespace mln;

  image2d<int> ima(5,5);

  iota(ima, 0);
  auto d = where(ima < 10);

  mln_foreach(const point2d& p, d)
    {
      BOOST_CHECK(d.has(p));
      BOOST_CHECK(ima(p) < 10);
    }

  mln_foreach(const auto& px, ima.pixels())
    {
      if (px.val() >= 10)
        BOOST_CHECK(not d.has(px.point()));
      else
        BOOST_CHECK(d.has(px.point()));
    }
}

BOOST_AUTO_TEST_CASE(Where_vfunction)
{
  using namespace mln;

  image2d<int> ima(5,5);

  iota(ima, 0);
  auto d = where(ima, [] (long x) { return x < 10; });

  mln_foreach(const point2d& p, d)
    {
      BOOST_CHECK(d.has(p));
      BOOST_CHECK(ima(p) < 10);
    }

  mln_foreach(const auto& px, ima.pixels())
    {
      if (px.val() >= 10)
        BOOST_CHECK(not d.has(px.point()));
      else
        BOOST_CHECK(d.has(px.point()));
    }
}



BOOST_AUTO_TEST_CASE(Where_pixfunction)
{
  using namespace mln;

  typedef  image2d<int> I;
  image2d<int> ima(5,5);

  iota(ima, 0);
  auto d = where(ima, [] (const mln_pixel(const I)& x) { return x.val() < 10; });

  mln_foreach(const point2d& p, d)
    {
      BOOST_CHECK(d.has(p));
      BOOST_CHECK(ima(p) < 10);
    }

  mln_foreach(const auto& px, ima.pixels())
    {
      if (px.val() >= 10)
        BOOST_CHECK(not d.has(px.point()));
      else
        BOOST_CHECK(d.has(px.point()));
    }
}
