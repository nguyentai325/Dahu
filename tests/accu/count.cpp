#include <mln/accu/accumulators/count.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Sum)
{
  using namespace mln::accu;

 accumulators::count<unsigned> acc;
 accumulators::count<unsigned> acc2;

 acc.take(12);
 acc.take("blabla");

 BOOST_CHECK_EQUAL( extractor::count(acc), 2);

 acc.untake(13);

 BOOST_CHECK_EQUAL( extractor::count(acc), 1);

 acc2.take(69);
 acc.take(69);
 acc2.take(acc);

 BOOST_CHECK_EQUAL( extractor::count(acc2), 3);
}
