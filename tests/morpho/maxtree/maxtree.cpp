#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/range/algorithm/generate.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/morpho/maxtree/maxtree.hpp>
#include <mln/core/grays.hpp>
#include <mln/io/imprint.hpp>

#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>
#include <random>


BOOST_AUTO_TEST_CASE(Maxtree)
{
  using namespace mln;
  typedef UInt<8> V;
  typedef typename image2d<V>::size_type size_type;
  image2d<V> ima(500, 500);


  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> sampler(0, value_traits<V>::max());
  range::generate(ima.values(), [&sampler, &gen] () { return sampler(gen); }) ;

  //iota(ima, 0);

  auto ctree = morpho::maxtree_indexes(ima, c4);

}
