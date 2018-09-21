#include <mln/core/win2d.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>



BOOST_AUTO_TEST_CASE(win2d_general)
{
  using namespace mln;
  {
    rect2d win = make_rectangle2d(3,5);

    mln_iter(p, win(literal::zero));
    p.init();
    for (int i = -1; i <= 1; ++i)
      for (int j = -2; j <= 2; ++j, p.next())
	BOOST_CHECK_EQUAL(*p, point2d(i,j));
  }

}
