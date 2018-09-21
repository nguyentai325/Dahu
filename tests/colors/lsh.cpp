#include <mln/core/colors.hpp>
#include <mln/colors/lsh.hpp>


#define BOOST_TEST_MODULE Colors
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(LSH)
{
  using namespace mln;

  typedef vec3i V;
  {
    rgb8 v = {0,0,0};
    BOOST_CHECK_EQUAL((V) v, (V) lsh2rgb(rgb2lsh(v)));
  }

  {
    rgb8 v = {255,255,255};
    BOOST_CHECK_EQUAL( (V) v,  (V) lsh2rgb(rgb2lsh(v)));
  }

  {
    rgb8 v = {255,128,0};
    BOOST_CHECK_EQUAL( (V) v,  (V) lsh2rgb(rgb2lsh(v)));
  }

  {
    lsh8 v = {0,0,128}; // Hue not significant
    BOOST_CHECK_EQUAL( (int) v[0],  (int) rgb2lsh(lsh2rgb(v))[0]);
  }

  {
    lsh8 v = {255,0,128}; // Hue and sat not significant
    BOOST_CHECK_EQUAL( (int) v[0],  (int) rgb2lsh(lsh2rgb(v))[0]);
  }

  {
    lsh8 v = {128,128,251}; // Hue on 252 values
    BOOST_CHECK_EQUAL( (V) v,  (V) rgb2lsh(lsh2rgb(v)));
  }

}
