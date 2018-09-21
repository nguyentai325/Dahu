#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/transform/chamfer_distance_transform.hpp>

#define BOOST_TEST_MODULE Transform
#include <tests/test.hpp>

BOOST_AUTO_TEST_SUITE(chamfer_distance_transform)

BOOST_AUTO_TEST_CASE(chamfer_distance_transform_1)
{
  using namespace mln;

  image2d<bool> f = { {0,0,0,0,0},
                      {0,1,1,1,0},
                      {0,1,1,1,0},
                      {0,1,1,1,0},
                      {0,0,0,0,0} };

  image2d<int> ref = { {0,0,0,0,0},
                       {0,1,1,1,0},
                       {0,1,2,1,0},
                       {0,1,1,1,0},
                       {0,0,0,0,0} };


  auto res = transform::chamfer_distance_transform(f, c4);

  MLN_CHECK_IMEQUAL(res, ref);
}

BOOST_AUTO_TEST_CASE(chamfer_distance_transform_2)
{
  using namespace mln;

  image2d<bool> f = { {1,0,0,0,0,1},
                      {0,1,1,1,1,0},
                      {0,1,1,1,1,0},
                      {0,1,1,1,1,0},
                      {1,0,0,0,0,1} };

  image2d<int> ref = { {1,0,0,0,0,1},
                       {0,1,1,1,1,0},
                       {0,1,2,2,1,0},
                       {0,1,1,1,1,0},
                       {1,0,0,0,0,1}};


  auto res = transform::chamfer_distance_transform(f, c4);

  MLN_CHECK_IMEQUAL(res, ref);
}


// Input and output are the same image
BOOST_AUTO_TEST_CASE(chamfer_distance_transform_3)
{
  using namespace mln;

  image2d<int> f = { {1,0,0,0,0,1},
                     {0,1,1,1,1,0},
                     {0,1,1,1,1,0},
                     {0,1,1,1,1,0},
                     {1,0,0,0,0,1} };

  image2d<int> ref = { {1,0,0,0,0,1},
                       {0,1,1,1,1,0},
                       {0,1,2,2,1,0},
                       {0,1,1,1,1,0},
                       {1,0,0,0,0,1}};


  transform::chamfer_distance_transform(f, c4, f);

  MLN_CHECK_IMEQUAL(f, ref);
}


BOOST_AUTO_TEST_SUITE_END()
