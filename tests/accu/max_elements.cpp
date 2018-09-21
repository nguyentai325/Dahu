#include <mln/accu/accumulators/max_elements.hpp>
#include <mln/core/colors.hpp>

#define BOOST_TEST_MODULE Accu
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(Max_elements)
{
  using namespace mln::accu;

 accumulators::maximal_elements<unsigned> acc;

 acc.take(5);
 acc.take(15);
 acc.take(2);

 {
   auto res = extractor::maximal_elements(acc);
   BOOST_CHECK_EQUAL(res.size(), 1u);
   BOOST_CHECK_EQUAL(res[0], 15u);
 }

 acc.take(15);
 acc.take(2);

 {
   auto res = extractor::maximal_elements(acc);
   BOOST_CHECK_EQUAL(res.size(), 1u);
   BOOST_CHECK_EQUAL(res[0], 15u);
 }

}

BOOST_AUTO_TEST_CASE(Max_elements_vec)
{
  using namespace mln;
  using namespace mln::accu;

  accumulators::maximal_elements<rgb8> acc;

  acc.take(rgb8 {(uint8)5, (uint8)5, (uint8)255});
  acc.take(rgb8 {(uint8)255, (uint8)5, (uint8)255});
  acc.take(rgb8 {(uint8)2, (uint8)2, (uint8)25});

  {
    auto res = extractor::maximal_elements(acc);
    BOOST_CHECK_EQUAL(res.size(), 1u);
    BOOST_CHECK_EQUAL(res[0], (rgb8 {(uint8)255, (uint8)5, (uint8)255}));
  }

  acc.take(rgb8 {(uint8)5, (uint8)8, (uint8)50});
  acc.take(rgb8 {(uint8)4, (uint8)1, (uint8)25});

  {
    auto res = extractor::maximal_elements(acc);
    BOOST_CHECK_EQUAL(res.size(), 2u);
    BOOST_CHECK_EQUAL(res[0], (rgb8 {(uint8)255, (uint8)5, (uint8)255}));
    BOOST_CHECK_EQUAL(res[1], (rgb8 {(uint8)5, (uint8)8, (uint8)50}));
  }

}
