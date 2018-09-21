#include <array>
#include <mln/core/range/zip.hpp>
#include <mln/core/forall.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>
#include <iostream>


BOOST_AUTO_TEST_CASE(ziprange)
{

  using namespace mln;

  std::array<int, 4> x {{2, 5, 15, 22}};
  std::array<int, 4> y {{-2, -5, -15, -22}};

  BOOST_CHECK_EQUAL(rng::size(x), 4);

  auto z = rng::zip(x, y);


  {
    int i = 0;
    mln_foreach(auto v, z) {
      BOOST_CHECK_EQUAL(std::get<0>(v), x[i]);
      BOOST_CHECK_EQUAL(std::get<1>(v), y[i++]);
    }
  }

  {
    int i = 0;
    mln_iter(pv, z);
    mln_forall(pv) {
      BOOST_CHECK_EQUAL(std::get<0>(*pv), x[i]);
      BOOST_CHECK_EQUAL(std::get<1>(*pv), y[i++]);
    }
  }
}
