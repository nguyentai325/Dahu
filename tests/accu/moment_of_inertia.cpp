#include <mln/core/point.hpp>
#include <mln/accu/accumulators/moment_of_inertia.hpp>

#define BOOST_TEST_MODULE Accu
#include <tests/test.hpp>


// Moment of inertia of a square
BOOST_AUTO_TEST_CASE(Moment_of_inertia)
{
  using namespace mln;



  auto sum_of_square = [](float x) { return x*(x+1)*(2*x+1)/6; };

  double res1;
  {
    accu::accumulators::moment_of_inertia<point2d, vec2u> acc;
    acc.init();
    int n = 10;
    for (short i = 0; i <= n; ++i)
      for (short j = 0; j <= n; ++j)
        acc.take(point2d{i,j});
    res1 = acc.to_result();

    double tmp = 2 * (n+1) * sum_of_square(n/2);
    double ref = (tmp + tmp) / sqr(sqr(n+1));
    BOOST_CHECK_EQUAL(res1, ref);
  }

  // We are supposed to be scale invariant.
  {
    accu::accumulators::moment_of_inertia<point2d, vec<double,2> > acc;
    acc.init();
    int n = 10000;
    for (short i = 0; i <= n; ++i)
      for (short j = 0; j <= n; ++j)
        acc.take(point2d{i,j});

    double res2 = acc.to_result();
    double tmp = 2 * (n+1) * sum_of_square(n/2);
    double ref = (tmp + tmp) / sqr(sqr(double(n+1)));
    BOOST_CHECK_CLOSE(res2, ref, 1);
    BOOST_CHECK_CLOSE(res1, res2, 1);
  }
}


