#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/range/algorithm/generate.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/io/imread.hpp>
#include <mln/morpho/alphatree/alphatree.hpp>
#include <mln/morpho/datastruct/checktree.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>
#include <mln/core/grays.hpp>


#define BOOST_TEST_MODULE Morpho
#include <tests/test.hpp>
#include <random>

/*
BOOST_AUTO_TEST_CASE(AlphaTree)
{
  using namespace mln;
  typedef uint8 V;
  typedef typename image2d<V>::size_type size_type;
  image2d<V> ima(100, 100);


  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> sampler(0, value_traits<V>::max());
  range::generate(ima.values(), [&sampler, &gen] () { return sampler(gen); }) ;

  //iota(ima, 0);

  typedef morpho::component_tree< unsigned, image2d<unsigned> > tree_t;
  tree_t                                tree;
  property_map<tree_t, int>             vmap;

  std::tie(tree, vmap) = morpho::alphatree_indexes(ima, c4);

  morpho::internal::checktree(tree);
}
*/

BOOST_AUTO_TEST_CASE(AlphaTree_2)
{
  using namespace mln;
  typedef uint8 V;
  typedef typename image2d<V>::size_type size_type;
  image2d<V> ima;

  io::imread(MLN_IMG_PATH "fly.pgm", ima);


  typedef morpho::component_tree< unsigned, image2d<unsigned> > tree_t;
  tree_t                                tree;
  property_map<tree_t, int>             vmap;

  std::tie(tree, vmap) = morpho::alphatree_indexes(ima, c4);

  morpho::internal::checktree(tree);
  std::cout << "Size: " << tree.realsize() << std::endl;

  {
    std::ofstream out("out.dot");
    morpho::write_graphviz(out, tree, vmap);
  }

}
