#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/vec/vec_io.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>

#include <mln/accu/accumulators/accu_if.hpp>
#include <mln/accu/accumulators/count.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <fstream>

#include <boost/format.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/transpose_graph.hpp>
//#include <boost/random/mersenne_twister.hpp>

# include <mln/morpho/component_tree/pattern_spectra.hpp>
# include <mln/morpho/component_tree/compute_depth.hpp>

#include "compute_g2.hpp"
#include "remove_parent_relation.hpp"
#include "reconstruct.hpp"
#include "routines.hpp"

namespace mln
{

  template <typename Graph, typename V, class GraphValueMap>
  void write_graphviz(const std::string& filename,
                      const tree_t& t1, const tree_t& t2, const tree_t& t3,
                      const property_map<tree_t, typename Graph::vertex_descriptor>& t1link,
                      const property_map<tree_t, typename Graph::vertex_descriptor>& t2link,
                      const property_map<tree_t, typename Graph::vertex_descriptor>& t3link,
                      const property_map<tree_t, V>& val1,
                      const property_map<tree_t, V>& val2,
                      const property_map<tree_t, V>& val3,
                      const Graph& graph,
                      const GraphValueMap& vmap)
  {
    using mln::io::format;

    /// Write the three trees
    {
      auto prop = [&](const tree_t::node_type& x) {
        return (boost::format("%i\\nval:%i") % t1link[x] % (int)val1[x]).str();
      };
      std::ofstream outs(filename + "-t1.dot");
      morpho::write_graphviz(outs, t1, make_functional_property_map(prop));
      outs.close();
    }
    {
      auto prop = [&](const tree_t::node_type& x) {
        return (boost::format("%i\\nval:%i") % t2link[x] % (int)val2[x]).str();
      };
      std::ofstream outs(filename + "-t2.dot");
      morpho::write_graphviz(outs, t2, make_functional_property_map(prop));
      outs.close();
    }
    {
      auto prop = [&](const tree_t::node_type& x) {
        return (boost::format("%i\\nval:%i") % t3link[x] % (int)val3[x]).str();
      };
      std::ofstream outs(filename + "-t3.dot");
      morpho::write_graphviz(outs, t3, make_functional_property_map(prop));
      outs.close();
    }
    /// Write the graph
    {
      auto vpwriter = [&vmap] (std::ostream& os, const typename Graph::vertex_descriptor& v) {
        os << "[label=\"" << v << "\\n"; format(os, vmap[v]) << "\"]";
      };

      std::ofstream outs(filename + ".dot");
      boost::write_graphviz(outs, graph, vpwriter);
      outs.close();
    }
  }
}


template <class WeightPropertyMap>
void
make_minimum_spanning_tree(Graph& g, const WeightPropertyMap& weights)
{
  mln_entering("Minimum Spanning tree computation");
  Graph tg;
  boost::transpose_graph(g, tg);

  auto idxpmap = boost::get(boost::vertex_index, tg);
  boost::vector_property_map<Graph::vertex_descriptor, decltype(idxpmap)> parent(boost::num_vertices(tg), idxpmap);

  boost::prim_minimum_spanning_tree(tg, parent,
                                    (boost::root_vertex(boost::vertex(0, tg)).
                                     weight_map(weights)));

  std::cout << "Before: num edges " << boost::num_edges(g) << std::endl;
  // Update g2
  boost::remove_edge_if([&tg, &parent](const Graph::edge_descriptor& e) {
      return parent[boost::source(e, tg)] != boost::target(e, tg); },
    g);

  std::cout << "aFTER: num edges " << boost::num_edges(g) << std::endl;
  mln_exiting();
}


template <class WeightPropertyMap, class Compare = std::less<typename boost::property_traits<WeightPropertyMap>::value_type>  >
void
make_shortest_path_dag(Graph& g, const WeightPropertyMap& weights, Compare cmp = Compare())
{
  mln_entering("Minimum Spanning tree computation");
  Graph tg;
  boost::transpose_graph(g, tg);

  auto idxpmap = boost::get(boost::vertex_index, tg);
  boost::vector_property_map<Graph::vertex_descriptor, decltype(idxpmap)> parent(boost::num_vertices(tg), idxpmap);

  typedef typename boost::property_traits<WeightPropertyMap>::value_type distance_type;

  boost::dag_shortest_paths(tg, boost::vertex(0, tg),
                            boost::predecessor_map(parent).
                            weight_map(weights).
                            distance_compare(cmp).
                            distance_inf(mln::value_traits<distance_type>::sup()).
                            distance_zero(distance_type(0))
                            );


  std::cout << "Before: num edges " << boost::num_edges(g) << std::endl;
  // Update g2
  boost::remove_edge_if([&tg, &parent](const Graph::edge_descriptor& e) {
      return parent[boost::source(e, tg)] != boost::target(e, tg); },
    g);

  std::cout << "aFTER: num edges " << boost::num_edges(g) << std::endl;
  mln_exiting();
}


template <class WeightPropertyMap, class Function>
float
compute_graph_energy(const Graph& g, const WeightPropertyMap& w, Function f, float init)
{
  Graph::edge_iterator it, end;
  std::tie(it,end) = boost::edges(g);
  float res = std::accumulate(it, end,
                              init,
                              [&w, f](float v, Graph::edge_descriptor e) { return f(v, w[e]); });

  return res;
}


boost::vector_property_map<unsigned>
compute_graph_depth(const Graph& g)
{
  boost::vector_property_map<unsigned> depth(boost::num_vertices(g));
  //auto viz_ = boost::record_distances(depth, boost::on_tree_edge ());
  //auto viz = boost::make_bfs_visitor(viz_);

  auto one = [](mln::dontcare_t) -> int{ return 1; };
  auto w = boost::make_function_property_map<Graph::edge_descriptor, int, decltype(one)>(one);

  Graph::vertex_descriptor root = boost::vertex(0, g);
  depth[root] = 0;

  Graph gT; 
  boost::transpose_graph(g, gT);
  boost::dag_shortest_paths(gT, root, boost::weight_map(w)
                            .distance_map(depth)
                            .distance_compare(std::greater<int> ())
                            .distance_inf(-1)
                            .distance_zero(0)
                            );
  //boost::breadth_first_search(gT, root, boost::visitor(viz));

  return depth;
}

template <class ValueMap, class BinaryFunction, class ValueType>
mln::image2d<ValueType>
write_vmap_to_image(const Graph& g, const mln::image2d<mln::vec3u>& pixmap,
                    const ValueMap& vmap, BinaryFunction op, ValueType init)
{
  (void) g;

  using namespace mln;
  image2d<ValueType> out;
  resize(out, pixmap);

  mln_pixter(px, out);
  mln_forall(px)
    {
      vec3u v = pixmap[px->index()];
      ValueType w = init;
      for (int k = 0; k < 3; ++k)
        w = op(w, vmap[v[k]]);
      px->val() = w;
    }

  return out;
}


// template <class NodePropertyMap, class EdgePropertyMap>
// mln::image2d<unsigned>
// compute_pattern_spectra(const Graph& g, const NodePropertyMap& pmap1, const EdgePropertyMap& pmap2)
// {
//   float amin, amax;
//   float wmin, wmax;
//   {
//     accu::accumulators::minmax<float> accu;
//     BOOST_FOREACH(Graph::vertex_descriptor v, boost::vertices(g)) {
//       accu.take(pmap1[v]);
//     }
//     std::tie(amin, amax) = accu.to_result();

//     accu.init();
//     BOOST_FOREACH(Graph::edge_descriptor e, boost::edges(g)) {
//       accu.take(pmap2[e]);
//     }
//     std::tie(wmin, wmax) = accu.to_result();
//   }

//   int nr = (int) log(amax) + 1;
//   int nc = (int) wmax + 1;

// }



int main(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm grain"
		<< std::endl;
      std::exit(1);
    }

  using namespace mln;
  typedef rgb8 value_t;

  image2d<rgb8> ima;
  io::imread(argv[1], ima);


  image2d<rgb8> f = addborder(ima, lexicographicalorder_less<rgb8>());
  image2d<rgb8> F = immerse_k1(f);

  typedef uint8 V;
  image2d<V> r = transform(f, [](rgb8 x) -> V { return x[0]; });
  image2d<V> g = transform(f, [](rgb8 x) -> V { return x[1]; });
  image2d<V> b = transform(f, [](rgb8 x) -> V { return x[2]; });

  /// Compute the marginal ToS
  tree_t trees[NTREE];
  for (int i = 0; i < NTREE; ++i)
    trees[i] = morpho::cToS( imtransform(f, [i](value_t x) -> V { return x[i]; }), c4);

  tree_t& t1 = trees[0];
  tree_t& t2 = trees[1];
  tree_t& t3 = trees[2];

  // 0. Debug pw min/max depth
  {
    auto d1 = morpho::compute_depth(t1);
    auto d2 = morpho::compute_depth(t2);
    auto d3 = morpho::compute_depth(t3);


    image2d<uint8> a,b,c;
    resize(a, F);
    resize(b, F);
    resize(c, F);
    reconstruction(t1, d1, a);
    reconstruction(t2, d2, b);
    reconstruction(t3, d3, c);
    auto x = imin(imin(a,b),c);
    auto y = imax(imax(a,b),c);
    io::imsave(eval(y), "depthmax.tiff");
    io::imsave(eval(x), "depthmin.tiff");
  }


  // {
  //   auto depth = morpho::compute_depth(t1);
  //   auto area = morpho::accumulate(t1, accu::features::count<>());
  //   image2d<float> out(200,500);
  //   morpho::pattern_spectra(t1, depth, area, out, false, true);
  //   io::imsave(transform(out, [](float x) { return std::log(1+x); }), "spectra.tiff");
  // }


  /// Compute the graph
  Graph g2;
  std::array< property_map<tree_t, Graph::vertex_descriptor>, NTREE> tlink;
  std::tie(g2, tlink) = compute_g2(trees);
  property_map<tree_t, Graph::vertex_descriptor>& t1link = tlink[0];
  property_map<tree_t, Graph::vertex_descriptor>& t2link = tlink[1];
  property_map<tree_t, Graph::vertex_descriptor>& t3link = tlink[2];


  // The marginal values of the ToS
  auto val1 = make_attribute_map_from_image(t1, immerse_k1(r));
  auto val2 = make_attribute_map_from_image(t2, immerse_k1(g));
  auto val3 = make_attribute_map_from_image(t3, immerse_k1(b));

  // report a per node color
  auto colors = compute_graph_node_colors(g2, val1, val2, val3);

  // compute the map pixel -> graph node
  image2d<vec3u> gmap;
  resize(gmap, F);
  compute_graph_map(g2, t1, t2, t3, t1link, t2link, t3link, gmap);

  // Save depth images
  {
    auto gdepth = compute_graph_depth(g2);
    image2d<unsigned> mindepth = write_vmap_to_image(g2, gmap, gdepth,
                                                     (const unsigned& (*) (const unsigned&, const unsigned&)) std::min<unsigned>,
                                                     value_traits<unsigned>::max());
    image2d<unsigned> maxdepth = write_vmap_to_image(g2, gmap, gdepth,
                                                     (const unsigned& (*) (const unsigned&, const unsigned&)) std::max<unsigned>,
                                                     0u);
    io::imsave(eval(imcast<uint16>(mindepth)), "before-mindepth.tiff");
    io::imsave(eval(imcast<uint16>(maxdepth)), "before-maxdepth.tiff");
  }

  auto glink = boost::get(&graph_content::tlinks, g2);

  /// Make a tree from the graph
  {
    // Put the weight on the edges: red edge -> -1 / GB edge -> 0
    // The weight represente the (inverse) depth of nearest red node
    typedef float R;

    auto edgefun = [&colors, &g2, &glink](const Graph::edge_descriptor& e) -> R {
      //return l2norm(colors[boost::source(e, g2)] - colors[boost::target(e, g2)]);
      return -1;
      // return vec3i{- (int) (glink[boost::source(e, g2)][0].id() != tree_t::npos()),
      //              - (int) (glink[boost::source(e, g2)][1].id() != tree_t::npos()),
      //              - (int) (glink[boost::source(e, g2)][2].id() != tree_t::npos()) };
    };
    boost::function_property_map<decltype(edgefun), Graph::edge_descriptor, R>  weights(edgefun);
    auto mycmp = [] (vec3i x, vec3i y) { return x[0] < y[0]; };

    //write_graphviz("before", t1, t2, t3, t1link, t2link, t3link, val1, val2, val3, g2, colors);
    float e1 = compute_graph_energy(g2, weights, (const float& (*)(const float&, const float&)) (std::max<float>), 0.0);
    std::cout << "Energy before:" << e1  << std::endl;

    //make_minimum_spanning_tree(g2, weights); //, mycmp);
    make_shortest_path_dag(g2, weights);
    update_parent_relation(g2, t1, t2, t3, t1link, t2link, t3link);

    float e2 = compute_graph_energy(g2, weights, (const float& (*)(const float&, const float&)) (std::max<float>), 0.0);
    std::cout << "Energy after:" << e2  << std::endl;
    //write_graphviz("after", t1, t2, t3, t1link, t2link, t3link, val1, val2, val3, g2, colors);
  }

  // Save depth images
  {
    auto gdepth = compute_graph_depth(g2);
    image2d<unsigned> mindepth = write_vmap_to_image(g2, gmap, gdepth,
                                                     (const unsigned& (*) (const unsigned&, const unsigned&)) std::min<unsigned>,
                                                     value_traits<unsigned>::max());
    image2d<unsigned> maxdepth = write_vmap_to_image(g2, gmap, gdepth,
                                                     (const unsigned& (*) (const unsigned&, const unsigned&)) std::max<unsigned>,
                                                     0u);

    io::imsave(eval(imcast<uint16>(mindepth)), "after-mindepth.tiff");
    io::imsave(eval(imcast<uint16>(maxdepth)), "after-maxdepth.tiff");
  }

  // compute_graph_map(g2, t1, t2, t3, t1link, t2link, t3link, gmap);
  // io::imsave(eval(imcast<rgb8>(gmap)), "after-gmap.tiff");

  // Filter graph
  accu::accumulators::accu_if<accu::accumulators::count<>, K1::is_face_2_t, point2d> acc;
  auto a1 = morpho::paccumulate(t1, F, acc);
  auto a2 = morpho::paccumulate(t2, F, acc);
  auto a3 = morpho::paccumulate(t3, F, acc);

  // 4. reconstruction
  t1._reorder_pset();
  t2._reorder_pset();
  t3._reorder_pset();

  image2d<rgb8> out, tmp;
  out.resize(box2d{f.domain().pmin, f.domain().pmax * 2 - 1});
  resize(tmp, ima);

  for (int i = 2; i < argc; ++i)
    {
      unsigned grain = std::atoi(argv[i]);
      auto gpred = [&,grain] (Graph::vertex_descriptor v) {
        return
        (glink[v][0] != t1.nend() and a1[glink[v][0]] > grain) or
        (glink[v][1] != t2.nend() and a2[glink[v][1]] > grain) or
        (glink[v][2] != t3.nend() and a3[glink[v][2]] > grain);
      };

      boost::filtered_graph<Graph, boost::keep_all, std::function<bool(Graph::vertex_descriptor)> > fg(g2, boost::keep_all(), gpred);


      // {
      //   auto pred = make_functional_property_map([grain, &a2](const tree_t::node_type& node) {
      //       return a2[node] >= grain; });
      //   filter_direct_and_reconstruct(t2, pred, val2, green(out));
      //   copy(green(out) | sbox2d{out.domain().pmin + 2, out.domain().pmax - 2, point2d{2,2}}, tmp);
      //   io::imsave(tmp, (boost::format("green-%06i.tiff") % grain).str().c_str());
      // }

      //reconstruct_marginal(t1, t2, t3, val1, val2, val3, fg, out);
      out = reconstruct_from_graph(g2, fg, t1, t2, t3, val1, val2, val3, gmap);
      copy(out | sbox2d{out.domain().pmin + 2, out.domain().pmax - 2, point2d{2,2}}, tmp);
      io::imsave(tmp, (boost::format("out-%06i.tiff") % grain).str().c_str());
    }

}
