#include <mln/accu/accumulators/mean.hpp>
#include <mln/core/colors.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Mean)
{
  using namespace mln::accu;

 accumulators::mean<unsigned> acc;
 accumulators::mean<unsigned> acc2;

 acc.take(5);
 acc.take(15);

 BOOST_CHECK_EQUAL( extractor::sum(acc), 20);
 BOOST_CHECK_EQUAL( extractor::count(acc), 2);
 BOOST_CHECK_EQUAL( extractor::mean(acc), 10);

 acc.untake(10);

 BOOST_CHECK_EQUAL( extractor::sum(acc), 10);
 BOOST_CHECK_EQUAL( extractor::count(acc), 1);
 BOOST_CHECK_EQUAL( extractor::mean(acc), 10);


 acc2.take(20);
 acc2.take(40);
 acc2.take(acc);

 BOOST_CHECK_EQUAL( extractor::sum(acc2), 70);
 BOOST_CHECK_EQUAL( extractor::count(acc2), 3);
 BOOST_CHECK_EQUAL( extractor::mean(acc2), 23);
}

BOOST_AUTO_TEST_CASE(Mean_vec)
{
  using namespace mln;
  using namespace mln::accu;

 accumulators::mean<rgb8> acc;
 accumulators::mean<rgb8> acc2;

 acc.take(rgb8 {(uint8)5, (uint8)5, (uint8)255});
 acc.take(rgb8 {(uint8)255, (uint8)5, (uint8)255});

 BOOST_CHECK_EQUAL( extractor::sum(acc), (rgb<int> {260,10,510})  );
 BOOST_CHECK_EQUAL( extractor::count(acc), 2);
 BOOST_CHECK_EQUAL( extractor::mean(acc), (rgb8 {(uint8)130, (uint8)5, (uint8)255}) );

 acc.untake(rgb8{(uint8)10,(uint8)10, (uint8)10} );

 BOOST_CHECK_EQUAL( extractor::sum(acc), (rgb<int> {250,0,500}) );
 BOOST_CHECK_EQUAL( extractor::count(acc), 1);
 BOOST_CHECK_EQUAL( extractor::mean(acc), (rgb<int> {250,0,500}) );
}
