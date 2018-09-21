#include <array>
#include <iostream>


#include <mln/morpho/datastruct/component_tree.hpp>
#include <mln/core/foreach.hpp>
#include <mln/core/forall.hpp>


#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(component_tree_structure)

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


BOOST_AUTO_TEST_CASE(component_tree_forward_iteration)
{
  tree_t tree = make_tree();

  unsigned n = 1;
  mln_iter(p, tree.pset());
  p.init();

  mln_foreach(auto x, tree.nodes())
    {
      BOOST_CHECK(x.id() == n++);
      std::cout << "== Node: " << x.id() << std::endl;
      std::cout << "{ first child: " << x.get_first_child_id()
		<< " next sibling: " << x.get_next_sibling_id()
		<< " parent: " << x.get_parent_id() << " }" << std::endl;
      std::cout << "Children: ";
      mln_foreach(auto y, x.children())
	std::cout << y.id() << ",";
      std::cout << std::endl;
      std::cout << "Proper pset: ";
      mln_foreach(auto y, x.proper_pset())
	{
	  BOOST_CHECK_EQUAL(y, *p);
	  std::cout << y << ",";
	  p.next();
	}
      std::cout << std::endl;
      std::cout << "Full pset: ";
      mln_foreach(auto y, x.pset())
	std::cout << y << ",";
      std::cout << std::endl;
    }
  BOOST_CHECK_EQUAL(n, 7);
}

BOOST_AUTO_TEST_CASE(component_tree_backward_iteration)
{
  tree_t tree = make_tree();

  unsigned n = 6;
  auto p = tree.pset().riter();
  p.init();
  mln_reverse_foreach(auto x, tree.nodes())
    {
      BOOST_CHECK(x.id() == n--);
      std::cout << "== Node: " << x.id() << std::endl;
      std::cout << "Proper pset: ";
      mln_reverse_foreach(auto y, x.proper_pset())
	{
	  BOOST_CHECK_EQUAL(y, *p);
	  std::cout << y << ",";
	  p.next();
	}
      std::cout << std::endl;
      std::cout << "Full pset: ";
      mln_reverse_foreach(auto y, x.pset())
	std::cout << y << ",";
      std::cout << std::endl;
    }
  BOOST_CHECK_EQUAL(n, 0);
}

BOOST_AUTO_TEST_CASE(component_tree_subtree)
{
  tree_t tree = make_tree();
  tree = tree.get_subtree(2);

  unsigned n = 0;
  mln_iter(p, tree.pset());
  p.init();

  mln_foreach(auto x, tree.nodes())
    {
      std::cout << "== Node: " << x.id() << std::endl;
      std::cout << "{ first child: " << x.get_first_child_id()
		<< " next sibling: " << x.get_next_sibling_id()
		<< " parent: " << x.get_parent_id() << " }" << std::endl;
      std::cout << "Children: ";
      mln_foreach(auto y, x.children())
	std::cout << y.id() << ",";
      std::cout << std::endl;
      std::cout << "Proper pset: ";
      mln_foreach(auto y, x.proper_pset())
	{
	  BOOST_CHECK_EQUAL(y, *p);
	  std::cout << y << ",";
	  p.next();
	}
      std::cout << std::endl;
      std::cout << "Full pset: ";
      mln_foreach(auto y, x.pset())
	std::cout << y << ",";
      std::cout << std::endl;
      n++;
    }
  BOOST_CHECK_EQUAL(n, 3);
}


BOOST_AUTO_TEST_SUITE_END()


