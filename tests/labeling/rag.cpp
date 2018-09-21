#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/labeling/rag.hpp>

#define BOOST_TEST_MODULE Labeling
#include <boost/test/unit_test.hpp>
#include <mln/io/imprint.hpp>


BOOST_AUTO_TEST_CASE(rag_test_1)
{
  using namespace mln;

  image2d<bool> ima(5,5);
  ima.at(0,0) = true;
  ima.at(0,4) = true;
  ima.at(4,0) = true;
  ima.at(4,4) = true;

  typedef boost::property<labeling::vertex_label_t, uint8> Vproperty;
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vproperty> Graph;

  image2d<uint8> lbl;
  unsigned nlabel;
  Graph graph;
  std::tie(lbl, nlabel) = labeling::rag(ima, c4, graph, uint8 ());

  io::imprint(lbl);
  {
    boost::graph_traits<Graph>::edge_iterator e,f;
    std::tie(e,f) = boost::edges(graph);
    for (; e != f; ++e)
      {
	auto x = boost::source(*e, graph);
	auto y = boost::target(*e, graph);
	std::cout << (int) boost::get(labeling::vertex_label_t (), graph, x) << "<-->"
		  << (int) boost::get(labeling::vertex_label_t (), graph, y) << std::endl;
      }
  }
}


BOOST_AUTO_TEST_CASE(rag_test_2)
{
  using namespace mln;

  image2d<bool> ima(5,5);
  ima.at(0,0) = true;
  ima.at(0,4) = true;
  ima.at(4,0) = true;
  ima.at(4,4) = true;

  typedef boost::property<labeling::vertex_label_t, uint8> Vproperty;
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vproperty> Graph;

  image2d<uint8> lbl;
  unsigned nlabel;
  Graph graph;
  std::tie(lbl, nlabel) = labeling::rag(ima, c8, graph, uint8 ());

  io::imprint(lbl);
  {
    boost::graph_traits<Graph>::edge_iterator e,f;
    std::tie(e,f) = boost::edges(graph);
    for (; e != f; ++e)
      {
	auto x = boost::source(*e, graph);
	auto y = boost::target(*e, graph);
	std::cout << (int) boost::get(labeling::vertex_label_t (), graph, x) << "<-->"
		  << (int) boost::get(labeling::vertex_label_t (), graph, y) << std::endl;
      }
  }
}
