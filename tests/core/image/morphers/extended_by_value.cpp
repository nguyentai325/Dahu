#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/image/morphers/extended_by_value_image.hpp>

#include <mln/core/algorithm/iota.hpp>

#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>




BOOST_AUTO_TEST_SUITE(ValueExtension)

BOOST_AUTO_TEST_CASE(value_extension_lvalue)
{
  using namespace mln;

  box2d dom{{-1,-2}, {3,3}};
  image2d<int> ima(dom);

  iota(ima, 0);
  auto x = extend_by_value(ima, 69);

  {
    mln_pixter(p, q, ima, x);
    mln_forall(p,q)
      BOOST_CHECK_EQUAL(p->val(), q->val());
  }

  {
    mln_pixter(p, x);
    mln_iter(q, c8(p));
    mln_forall(p)
      mln_forall(q)
    {
      if (!x.domain().has(q->point()))
	BOOST_CHECK_EQUAL(q->val(), 69);
      else
        BOOST_CHECK_EQUAL(&q->val(), &ima(q->point()));
    }
  }
}


BOOST_AUTO_TEST_CASE(value_extension_rvalue)
{
  using namespace mln;

  box2d dom{{-1,-2}, {3,3}};
  image2d<int> ima(dom);

  iota(ima, 0);
  auto x = extend_by_value(ima, 69);

  {
    mln_pixter(p, q, ima, x);
    mln_forall(p,q)
      BOOST_CHECK_EQUAL(p->val(), q->val());
  }

  {
    mln_pixter(p, x);
    mln_iter(q, c8(p));
    mln_forall(p)
      mln_forall(q)
    {
      if (!x.domain().has(q->point()))
	BOOST_CHECK_EQUAL(q->val(), 69);
      else
        BOOST_CHECK_EQUAL(&q->val(), &ima(q->point()));
    }
  }
}



BOOST_AUTO_TEST_SUITE_END()
