#include <mln/core/domain/box.hpp>


#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Box)

BOOST_AUTO_TEST_CASE(box2d_general)
{
  using namespace mln;
  {
    box2d a ({0,0}, {0,0});
    BOOST_CHECK(a.empty());
  }

  {
    short minp = - (1 << 15);
    short maxp = (1 << 15) - 1;
    unsigned n = (1<<16)-1;
    box2d a ({minp,minp}, {maxp,maxp});
    BOOST_CHECK_EQUAL(a.size(), n*n);
  }

}


BOOST_AUTO_TEST_CASE(box2d_forward)
{
  using namespace mln;

  box2d b({2,3}, {6,8});

  point2d p;
  int i = 2, j = 3;
  mln_foreach(p, b) {
    BOOST_CHECK_EQUAL(p[0], i);
    BOOST_CHECK_EQUAL(p[1], j);
    BOOST_CHECK(b.has(p));
    if (++j == 8) { j = 3; ++i; }
  }

}

BOOST_AUTO_TEST_CASE(box2d_backward)
{
  using namespace mln;

  box2d b({2,3}, {6,8});

  point2d p;
  int i = 5, j = 7;
  mln_foreach(p, b.riter()) {
    BOOST_CHECK_EQUAL(p[0], i);
    BOOST_CHECK_EQUAL(p[1], j);
    BOOST_CHECK(b.has(p));
    if (--j < 3) { j = 7; --i; }
  }
}

BOOST_AUTO_TEST_CASE(strided_box2d_forward)
{
  using namespace mln;

  sbox2d b({2,3}, {7,16}, {2,3});
  std::cout << b << std::endl;

  auto p = b.iter(); p.init();
  for (int i = 2; i < 7; i+=2)
    for (int j = 3; j < 16; j+=3)
      {
	BOOST_CHECK(!p.finished());
	BOOST_CHECK_EQUAL((*p)[0], i);
	BOOST_CHECK_EQUAL((*p)[1], j);
	BOOST_CHECK(b.has(*p));
	p.next();
      }
  BOOST_CHECK(p.finished());
}

BOOST_AUTO_TEST_CASE(strided_box2d_backward)
{
  using namespace mln;

  sbox2d b({2,3}, {7,17}, {2,3});

  auto p = b.riter(); p.init();
  for (int i = 6; i >= 2; i-=2)
    for (int j = 15; j >= 3; j-=3)
      {
	BOOST_CHECK(!p.finished());
	BOOST_CHECK_EQUAL((*p)[0], i);
	BOOST_CHECK_EQUAL((*p)[1], j);
	BOOST_CHECK(b.has(*p));
	p.next();
      }
  BOOST_CHECK(p.finished());
}


BOOST_AUTO_TEST_SUITE_END()
