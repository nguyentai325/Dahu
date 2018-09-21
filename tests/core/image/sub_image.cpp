#define BOOST_TEST_MODULE Core
#include <tests/test.hpp>

#include <mln/core/image/image2d.hpp>
#include <mln/core/image/sub_image.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/algorithm/fill.hpp>




BOOST_AUTO_TEST_SUITE(SubImage)

BOOST_AUTO_TEST_CASE(sub_domain_with_box)
{
  using namespace mln;

  image2d<int> ima(5,5);

  {
    image2d<int> ref = {{0,1,2,3,4},
                        {5,42,42,42,9},
                        {10,42,42,42,14},
                        {15,42,42,42,19},
                        {20,21,22,23,24}};

    iota(ima, 0);
    fill(ima | box2d{{1,1}, {4,4}}, 42);
    MLN_CHECK_IMEQUAL(ima, ref);
  }

  {
    image2d<int> ref = {{0,1,2,3,4},
                        {5,42,42,8,9},
                        {10,11,12,13,14},
                        {15,16,17,18,19},
                        {20,21,22,23,24}};

    iota(ima, 0);
    fill(ima | box2d{{1,1}, {4,4}} | box2d{{1,1}, {2,3}}, 42);
    MLN_CHECK_IMEQUAL(ima, ref);
  }

  static_assert( std::is_same< decltype(ima | box2d()), image2d<int> >::value, "");
}





BOOST_AUTO_TEST_CASE(sub_domain)
{
  using namespace mln;

  image2d<int> ima(5,5);

  {
    image2d<int> ref = {{0,1,2,3,4},
                        {5,6,7,8,9},
                        {10,42,42,42,42},
                        {42,42,42,42,42},
                        {42,42,42,42,42}};
    iota(ima, 0);
    fill(ima | where(ima > 10), 42);
    MLN_CHECK_IMEQUAL(ima, ref);
  }


  {
    image2d<int> ref = {{0,1,2,3,4},
                        {5,6,7,8,9},
                        {10,11,12,13,14},
                        {15,16,17,18,19},
                        {20,42,42,42,42}};

    iota(ima, 0);
    fill((ima | where(ima > 10)) | where(ima > 20), 42);
    MLN_CHECK_IMEQUAL(ima, ref);
  }

  {
    image2d<int> ref = {{0,1,2,3,4},
                        {5,6,7,8,9},
                        {10,42,42,42,42},
                        {42,42,42,42,42},
                        {20,21,22,23,24}};
    iota(ima, 0);
    fill(ima | where(land(ima > 10, ima < 20)), 42);
    MLN_CHECK_IMEQUAL(ima, ref);
  }

}

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(sub_domain_failure, 1)

BOOST_AUTO_TEST_CASE(sub_domain_failure)
{
  using namespace mln;

  image2d<int> ima(5,5);
  iota(ima, 0);
  fill(ima | where(ima > 10), 42);
  // This should throw because
  // where(ima < 20) is not included in where(ima > 10)
  {
    iota(ima, 0);
    fill((ima | where(ima > 10)) | where(ima < 20), 44);
  }
}


/*
BOOST_AUTO_TEST_CASE(mask)
{
  using namespace mln;

  image2d<int> ima(5,5);
  iota(ima, 0);
  fill(ima | (ima > 12), 12);
  fill(ima | (ima < 5) , 5);


  static_assert( std::is_same< decltype(ima | box2d()), image2d<int> >::value, "");
  io::imprint(ima);
}
*/


BOOST_AUTO_TEST_SUITE_END()


