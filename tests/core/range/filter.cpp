#include <mln/core/domain/box.hpp>
#include <mln/core/forall.hpp>
#include <mln/core/range/filter.hpp>


#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>
#include <iostream>


BOOST_AUTO_TEST_CASE(range_filter)
{

  using namespace mln;

  box2d a ({0,0}, {6,11});
  auto x = rng::filter(a, [](point2d p) { return (p[0] % 2) == (p[1] % 2); }); // chess board

  unsigned sz = 0;
  mln_foreach(point2d p, x) {
    (void) p;
    ++sz;
  }

  BOOST_CHECK(x.has(point2d{0,0}));
  BOOST_CHECK(x.has(point2d{1,1}));
  BOOST_CHECK(not x.has(point2d{0,1}));
  BOOST_CHECK(not x.has(point2d{1,0}));
  BOOST_CHECK(sz == 33);
}
