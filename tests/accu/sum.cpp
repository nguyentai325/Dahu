#include <mln/accu/accumulators/sum.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Sum)
{
  using namespace mln::accu;

 accumulators::sum<unsigned> acc;
 accumulators::sum<unsigned> acc2;

 acc.take(12);
 acc.take(13);

 BOOST_CHECK_EQUAL( extractor::sum(acc), 25);

 acc.untake(10);

 BOOST_CHECK_EQUAL( extractor::sum(acc), 15);

 acc2.take(69);
 acc2.take(acc);

 BOOST_CHECK_EQUAL( extractor::sum(acc2), 84);
}
