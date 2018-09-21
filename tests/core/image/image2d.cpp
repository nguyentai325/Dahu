#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/extension/fill.hpp>
#include <mln/core/extension/mirror.hpp>
#include <numeric>

#define BOOST_TEST_MODULE Core
#include <tests/test.hpp>

BOOST_AUTO_TEST_SUITE(_image2d)


BOOST_AUTO_TEST_CASE(extension_fill)
{
  using namespace mln;

  image2d<int> ima(5,5);
  image2d<int> ref(5,5);


  iota(ima, 1);
  iota(ref, 1);

  mln::extension::fill(ima, 69);
  BOOST_CHECK_EQUAL(ima.at(-1,-1), 69);

  mln_pixter(p1, p2, ima, ref);
  mln_iter(q1, c4(p1));
  mln_iter(q2, c4(p2));

  int i = 1;
  mln_forall(p1, p2)
    {
      BOOST_CHECK_EQUAL(p1->val(), p2->val());
      BOOST_CHECK_EQUAL(p1->val(), i++);
      mln_forall(q1, q2)
	{
	  if (ima.domain().has(q1->point()))
	    BOOST_CHECK_EQUAL(q1->val(), q2->val());
	  else
	    BOOST_CHECK_EQUAL(q1->val(), 69);
	}
    }
}

BOOST_AUTO_TEST_CASE(extension_mirror)
{
  using namespace mln;

  image2d<int> ima(5,5);
  iota(ima, 1);

  mln::extension::mirror(ima);
  box2d dom = {{-3,-3}, {8,8}};

  // MIRROR f(p) = f(n - |p % 2n - n|
  {
    mln_foreach(point2d p, dom)
      {
        point2d q = p % 10;
        q[0] = q[0] < 0 ? q[0] + 10 : q[0];
        q[1] = q[1] < 0 ? q[1] + 10 : q[1]; // because c++ modulo rounds toward 0 (also for neg) 
        q[0] = q[0] < 5 ? q[0] : 9 - q[0];
        q[1] = q[1] < 5 ? q[1] : 9 - q[1];

        BOOST_CHECK(ima.at(p) == ima.at(q));
      }
  }

}

BOOST_AUTO_TEST_CASE(from_buffer)
{
  using namespace mln;
  int buffer[12];
  std::iota(buffer, buffer + 12, 1);

  image2d<int> ima = image2d<int>::from_buffer(buffer, box2d{{0,0},{3,4}});
  image2d<int> ref(3,4);

  iota(ref, 1);

  MLN_CHECK_IMEQUAL(ima, ref);
}


BOOST_AUTO_TEST_SUITE_END()
