#include <mln/graphcut/graphcut.hh>
#include <mln/core/image/image2d.hpp>
#include <mln/core/image/sub_image.hpp>

#include <mln/core/neighb2d.hpp>
#include <mln/core/foreach.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/io/imprint.hpp>

#define BOOST_TEST_MODULE GraphCut
#include <boost/test/unit_test.hpp>



BOOST_AUTO_TEST_CASE(graph_iteration)
{
  using namespace mln;

  box2d b({-1,-2}, {2,3});

  typedef mln::graphcut::internal::graphcut_graph_t<char, float, c4_t> graph_t;
  graph_t g(b, c4);
  unsigned nv = 0, ne = 0;

  mln_foreach(point2d p, g.vertices()) {
    g.vertex(p) = 'a' + (nv++ % 26);
  }

  mln_foreach(const auto& p, g.edges())
    g.edge(p) = ne++;

  {
    std::cout << "Vertex list: ";
    mln_foreach(point2d p, g.vertices())
      std::cout << "(" << p << ":" << g.vertex(p) << "), ";
    std::cout << std::endl;
  }

  typedef graph_t::edge_type edge_t;

  auto printe = [] (const edge_t& x) { std::cout << x.type << '|' << x.p; };
  {
    std::cout << "Edges list: ";
    mln_foreach(auto p, g.edges()) {
      std::cout << "(";
      printe(p);
      std::cout << ":" << g.edge(p) << "), ";
    }
    std::cout << std::endl;
  }


  int n = b.shape()[0];
  int m = b.shape()[1];
  BOOST_CHECK_EQUAL(nv, 2+n*m);
  BOOST_CHECK_EQUAL(ne, 2*n*m + n*(m-1) + (n-1)*(2*m-1));
}

BOOST_AUTO_TEST_CASE(graph_nbh)
{
  using namespace mln;

  box2d b({-1,-2}, {2,3});

  typedef mln::graphcut::internal::graphcut_graph_t<char, float, c4_t> graph_t;
  graph_t g(b, c4);
  unsigned nv = 0, ne = 0;

  mln_foreach(point2d p, g.vertices())
    g.vertex(p) = 'a' + (nv++ % 26);

  mln_foreach(const auto& e, g.edges())
    g.edge(e) = ne++;


  {
    std::cout << "Vertex list: (adj vertex)" << std::endl;
    mln_foreach(point2d p, g.vertices())
      {
	std::cout << "(" << p << ":" << g.vertex(p) << "): ";
	mln_foreach(point2d q, g.adjacent_vertices(p))
	  std::cout << "(" << q << ":" << g.vertex(q) << "), ";
	std::cout << std::endl;
      }
    std::cout << std::endl;
  }

  typedef graph_t::edge_type edge_t;
  auto printe = [] (const edge_t& x) { std::cout << x.type << '|' << x.p; };
  {
    std::cout << "Vertex list (adj edges): " << std::endl;
    mln_foreach(point2d p, g.vertices())
      {
	std::cout << "(" << p << ":" << g.vertex(p) << "): ";
	mln_foreach(const edge_t& q, g.edges(p)) {
	  std::cout << "(";
	  printe(q);
	  std::cout << ":" << g.edge(q) << "), ";
	}
	std::cout << std::endl;
      }
    std::cout << std::endl;
  }

}

  BOOST_AUTO_TEST_CASE(graphcut_1)
  {
    using namespace mln;

    image2d<bool> ori(5,5);
    fill(ori, false);
    fill(ori | sbox2d({0,0},{5,5},{2,2}), true);


    image2d<bool> out(5,5);

    graphcut::graphcut(ori, out, c4,
		       [](bool x, bool y) { return (int) 2 * (x != y); },
		       [](bool x, bool y) { return (int) 1 * (x != y); });

    io::imprint(out);
  }


  BOOST_AUTO_TEST_CASE(graphcut_2)
  {
    using namespace mln;

    image2d<float> ori(5,5);
    fill(ori, 0.0);
    fill(ori | sbox2d({0,0},{5,5},{2,2}), 0.7);


    image2d<bool> out(5,5);

    graphcut::graphcut(ori, out, c4,
		       [](bool x, float v) { return  2.0 * std::abs(x - v); },
		       [](bool x, bool y) { return 0.2 * (x != y); });

    io::imprint(out);
  }

