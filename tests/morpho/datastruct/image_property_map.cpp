#include <mln/morpho/datastruct/component_tree.hpp>
#include <mln/morpho/datastruct/image_property_map.hpp>
#include <mln/core/foreach.hpp>
#include <mln/core/forall.hpp>
#include <mln/core/algorithm/iota.hpp>


#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(property_image)

using namespace mln;
typedef morpho::component_tree<char, std::array<char, 256> > tree_t;

tree_t make_tree()
{
  tree_t tree;

  tree_t::_data_t* data = tree._get_data();
  data->m_pset_ordered = true;

  // We encode this tree
  //                   E
  //      a,b,c
  //     /     ╲
  //   d,e      i,j
  //   / \       |
  //  f  g,h     k

  data->m_nodes = {
    { 0, 6, 1, 0, 11}, // sentinel
    { 0, 0, 2, 0, 0},
    { 1, 1, 3, 5, 3},
    { 2, 2, 4, 4, 5},
    { 2, 3, 5, 5, 6},
    { 1, 4, 6, 0, 8},
    { 5, 5, 0, 0, 10}
  };

  data->m_S = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k'};

  data->m_pmap['a'] = data->m_pmap['b'] = data->m_pmap['c'] = 1;
  data->m_pmap['d'] = data->m_pmap['e'] = 2;
  data->m_pmap['f'] = 3;
  data->m_pmap['g'] = data->m_pmap['h'] = 4;
  data->m_pmap['i'] = data->m_pmap['j'] = 5;
  data->m_pmap['k'] = 6;

  return tree.get_subtree(1);
}


BOOST_AUTO_TEST_CASE(property_image)
{
  tree_t tree = make_tree();
  property_map<tree_t, int> vmap(tree);

  auto ima = make_image(tree, vmap);
  iota(ima, 0);

  {
    int i = 0;
    mln_pixter(px, ima);
    mln_iter(nit, tree.nodes());
    mln_forall(px, nit)
    {
      BOOST_CHECK_EQUAL(px->val(), i++);
      BOOST_CHECK(px->point() == *nit);
    }
  }
}

BOOST_AUTO_TEST_CASE(property_image_nbh)
{
  tree_t tree = make_tree();
  property_map<tree_t, int> vmap(tree);

  auto ima = make_image(tree, vmap);
  iota(ima, 1);

  {
    std::vector<unsigned> voisins[7] = {
      {},
      {2,5},
      {1,3,4},
      {2},
      {2},
      {1,6},
      {5}
    };

    morpho::tree_neighb_t nbh;
    {
      mln_pixter(px, ima);
      mln_iter(nx, nbh(px));
      int i = 1;
      mln_forall(px)
      {
        unsigned j = 0;
        mln_forall(nx)
        {
          BOOST_CHECK_EQUAL(nx->val(), voisins[i][j]);
          BOOST_CHECK_EQUAL(nx->point().id(), voisins[i][j]);
          ++j;
        }
        ++i;
      }
    }

    {
      mln_piter(p, ima);
      mln_iter(n, nbh(p));
      int i = 1;
      mln_forall(p)
      {
        unsigned j = 0;
        mln_forall(n)
        {
          BOOST_CHECK_EQUAL(ima(*n), voisins[i][j]);
          BOOST_CHECK(n->id() == voisins[i][j]);
          ++j;
        }
        ++i;
      }
    }

  }
}



BOOST_AUTO_TEST_SUITE_END()
