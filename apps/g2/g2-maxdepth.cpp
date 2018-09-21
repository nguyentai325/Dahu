#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/dontcare.hpp>
#include <mln/core/vec_base.hpp>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>

#include <mln/accu/accumulators/count.hpp>
#include <mln/accu/accumulators/accu_if.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>

#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/foreach.hpp>

#include "compute_g2.hpp"
#include "routines.hpp"

namespace mln
{
  boost::vector_property_map<unsigned>
  compute_graph_depth(const MyGraph& g)
  {
    mln_entering("Compute graph depth");

    boost::vector_property_map<unsigned> depth(boost::num_vertices(g));

    auto one = [](mln::dontcare_t) -> int{ return 1; };
    auto w = boost::make_function_property_map<MyGraph::edge_descriptor, int, decltype(one)>(one);

    MyGraph::vertex_descriptor root = boost::vertex(0, g);
    depth[root] = 0;

    MyGraph gT;
    boost::transpose_graph(g, gT);
    boost::dag_shortest_paths(gT, root, boost::weight_map(w)
                              .distance_map(depth)
                              .distance_compare(std::greater<int> ())
                              .distance_inf(-1)
                              .distance_zero(0)
                            );


    mln_exiting();
    return depth;
  }

  
  boost::vector_property_map<unsigned>
  compute_graph_count(const MyGraph& g)
  {

    struct viz_t : public boost::default_dfs_visitor
    {
      viz_t(boost::vector_property_map<unsigned>& common)
        : m_common (common)
      {
      }

      void examine_edge(MyGraph::edge_descriptor e, const MyGraph& g)
      {
        MyGraph::vertex_descriptor s = boost::source(e, g);
        MyGraph::vertex_descriptor t = boost::target(e, g);
        m_common[t] += m_common[s];
      }

      boost::vector_property_map<unsigned>& m_common;
    };

    mln_entering("Compute graph count");

    // 1. Set for each node its redondancy (the number of tree it belongs to minus 1)
    boost::vector_property_map<unsigned> common(boost::num_vertices(g));
    {
      auto tlinks = boost::get(&my_graph_content::tlinks, g);
      BOOST_FOREACH(MyGraph::vertex_descriptor v, boost::vertices(g))
        {
          common[v] = 0;
          for (int i = 0; i < NTREE; ++i)
            if (tlinks[v][i].id() != tree_t::npos())
              common[v] += 1;
          common[v] -= 1;
        }
      MyGraph::vertex_descriptor root = 0;
      common[root] = 0;
    }

    // 2. Propagate i.e. for each node X, compute the number of redondancy in X↑
    MyGraph gT;
    boost::transpose_graph(g, gT);
    boost::depth_first_search(gT, boost::visitor(viz_t(common)));

    // 3. Set for each node, the number of node it is included in.
    boost::vector_property_map<unsigned>& count = common; // do not need another vector, inplace
    {
      auto depth = boost::get(&my_graph_content::depth, g);
      BOOST_FOREACH(MyGraph::vertex_descriptor v, boost::vertices(g))
        count[v] = sum(depth[v]) - common[v];
    }


    mln_exiting();
    return count;
  }

  /// \brief Compute for each point the number of shapes it belongs to.
  /// \[ ω(x) = { X ∈ S, x ∈ X } \]
  image2d<unsigned>
  compute_map_count(const MyGraph& g,
                    const tree_t* t,
                    const tlink_t* tlink)
  {
    // 1. Compute the weight of each node to remove the redondancy
    // e.g., if a node appears twice, its weight is 1/2
    // and we sum up the weights throurgh the paths
    // We first store in redondancy_map[X] the number of trees X belongs to.
    boost::vector_property_map<char> redondancy_map(boost::num_vertices(g));
    auto glinks = boost::get(&my_graph_content::tlinks, g);

    BOOST_FOREACH(MyGraph::vertex_descriptor v, boost::vertices(g))
      {
        redondancy_map[v] = 0;
        for (int i = 0; i < NTREE; ++i)
          redondancy_map[v] += (glinks[v][i].id() != t[i].npos());
      }

    property_map<tree_t, float> w[NTREE];
    for (int i = 0; i < NTREE; ++i)
      {
        w[i] = property_map<tree_t, float> (t[i]);

        w[i][t[i].get_root()] = 0;
        mln_foreach(auto x, t[i].nodes_without_root()) // downward
          {
            MyGraph::vertex_descriptor v = tlink[i][x];
            w[i][x] = w[i][x.parent()] + 1.0 / redondancy_map[v];
          }
      }

    // 2. We set the weights value point-wise
    image2d<unsigned> dmap;
    resize(dmap, t[0]._get_data()->m_pmap).init(0);

    mln_foreach(auto px, dmap.pixels())
      {
        float r = 0;
        for (int i = 0; i < NTREE; ++i)
          {
            tree_t::node_type x = t[i].get_node_at(px.index());
            r += w[i][x];
          }
        px.val() = std::lround(r);
      }

    return dmap;
  }


  template <typename V>
  boost::vector_property_map<unsigned>
  compute_graph_variation(const MyGraph& g,
                          const boost::vector_property_map< vec<V, NTREE> >& colors)
  {
    mln_entering("Compute graph variation");

    boost::vector_property_map<unsigned> variation(boost::num_vertices(g));


    auto e_val = [&g,&colors](MyGraph::edge_descriptor e) -> int {
      return linfnorm(colors[boost::source(e,g)] - colors[boost::target(e,g)]);
    };
    auto w = boost::make_function_property_map<MyGraph::edge_descriptor, int, decltype(e_val)>(e_val);

    MyGraph::vertex_descriptor root = boost::vertex(0, g);
    variation[root] = 0;

    MyGraph gT;
    boost::transpose_graph(g, gT);
    boost::dag_shortest_paths(gT, root, boost::weight_map(w)
                              .distance_map(variation)
                              .distance_compare(std::greater<int> ())
                              .distance_inf(-1)
                              .distance_zero(0)
                            );

    mln_exiting();
    return variation;
  }

  template <typename V>
  boost::vector_property_map<float>
  compute_graph_area(const MyGraph& g,
                     const image2d<V>& ima,
                     const tree_t* trees,
                     const tlink_t* tlinks)
  {
    mln_entering("Compute graph area");

    boost::vector_property_map<float> area(boost::num_vertices(g));

    typedef accu::accumulators::count<unsigned> ACCU;
    accu::accumulators::accu_if<ACCU, K1::is_face_2_t, point2d> accu;

    for (int i = 0; i < NTREE; ++i)
      {
        auto amap = morpho::paccumulate(trees[i], ima, accu);
        float maxv = amap[trees[i].get_root()];
        mln_foreach(tree_t::node_type x, trees[i].nodes())
          area[tlinks[i][x]] = std::log(maxv - amap[x] + 1.0) / std::log(maxv + 1.0);
      }

    mln_exiting();
    return area;
  }


  template <class ValueMap, class BinaryFunction, class ValueType>
  void
  write_vmap_to_image(const MyGraph& g, const tree_t* t, const tlink_t* tlink,
                      const ValueMap& vmap, BinaryFunction op, ValueType init,
                      image2d<ValueType>& out)
  {
    (void) g;

    mln_pixter(px, out);
    mln_forall(px)
    {
      ValueType w = init;
      for (int i = 0; i < NTREE; ++i)
        {
          tree_t::node_type tnode = t[i].get_node_at(px->index());
          MyGraph::vertex_descriptor gnode = tlink[i][tnode];
          w = op(w, vmap[gnode]);
        }
      px->val() = w;
    }

  }

  template <class ValueMap, class ColorMap, class Compare,
            class ValueType, class ColorType>
  void
  write_vmap_to_image_and_rec(const MyGraph& g, const tree_t* t, const tlink_t* tlink,
                              const ValueMap& vmap,
                              const ColorMap& cmap,
                              Compare cmp,
                              image2d<ValueType>& vout,
                              image2d<ColorType>& cout)
  {
    (void) g;

    mln_pixter(px, vout);
    mln_pixter(px2, cout);
    mln_forall(px, px2)
    {
      tree_t::node_type tnode = t[0].get_node_at(px->index());
      MyGraph::vertex_descriptor current = tlink[0][tnode];
      int from = 0;
      for (int i = 1; i < NTREE; ++i)
        {
          tree_t::node_type tnode = t[i].get_node_at(px->index());
          MyGraph::vertex_descriptor gnode = tlink[i][tnode];
          if (cmp(vmap[gnode], vmap[current])) {
            current = gnode;
            from = i;
          }
        }
      px->val() = vmap[current];
      px2->val() = cmap[current];
      //px2->val() = 0;
      //px2->val()[from] = 255;
    }
  }

}

void usage(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " attribute input.ppm output.tiff [out1.tiff out2.tiff out3.tiff]\n"
    "Compute the graph G₂ and outputs a 16-bit depth image.\n"
    "Attribute in {'variation', 'depth', 'count', 'area'}\n"
    "  Variation:    somme des variations le long du plus long chemin\n"
    "  Depth:        longueur du plus long chemin d'inclusion (nombre de ligne de niveau\n"
    "                qui s'incluent max traversées depuis chaque point)\n"
    "  Count:        nombre de lignes traversées (nombre de noeuds incluant)\n"
    "  Area:         aire\n"
    "If 'out' is supplied, it also compute the individual RGB depth map.\n"
    ;
  std::exit(1);
}


int main(int argc, char** argv)
{
  enum e_opt { OPT_DEPTH, OPT_COUNT, OPT_COUNT2, OPT_VARIATION, OPT_AREA };
  std::map<std::string, e_opt> str2opt = {
    {"depth", OPT_DEPTH},
    {"count", OPT_COUNT},
    {"count2", OPT_COUNT2},
    {"variation", OPT_VARIATION},
    {"area", OPT_AREA}
  };

  if (argc < 4 or str2opt.find(argv[1]) == str2opt.end())
    usage(argv);


  tbb::task_scheduler_init init;

  using namespace mln;

  e_opt attribute = str2opt[argv[1]];
  const char* inname = argv[2];
  const char* outname = argv[3];
  argv += 4;
  argc -= 4;


  image2d<value_t> ima;
  io::imread(inname, ima);

  image2d<value_t> f = addborder_marginal(ima);
  image2d<value_t> F = immerse_k1(f);

  typedef value_t::value_type V;

  /// Compute the marginal ToS
  bool ASK = false;
  point2d pmin{0,0};
  if (ASK) {
    std::cout << "Point a l'infini: " << std::endl;
    std::cin >> pmin[0] >> pmin[1];
  }

  tree_t t[NTREE];

  tbb::parallel_for(0, (int)NTREE, [&t,&f,pmin](int i){
      t[i] = morpho::cToS_pinf(imtransform(f, [i](value_t x) { return x[i]; }), c4, pmin);
    });

  /// Compute the graph
  MyGraph g2;
  std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
  std::tie(g2, tlink) = compute_g2<NTREE>(t);


  /// If we have selected the depth attribute
  if (attribute == OPT_DEPTH or attribute == OPT_COUNT2)
    {
      /// Compute depth
      boost::vector_property_map<unsigned> gdepth;
      if (attribute == OPT_DEPTH)
        gdepth = compute_graph_depth(g2);
      else if (attribute == OPT_COUNT2)
        gdepth = compute_graph_count(g2);

      /// Compute individual depth (for debugging)
      if (argc > 0)
        {
          //property_map<tree_t, unsigned> d;
          image2d<unsigned> d_;
          resize(d_, F);

          for (int i = 0; i < std::min<int>(argc, NTREE); ++i)
            {
              auto d = make_functional_property_map<tree_t::node_type>([&tlink, &gdepth,i](tree_t::node_type x) {
                  return gdepth[tlink[i][x]];
                });
              //d = morpho::compute_depth(t[i]);
              morpho::reconstruction(t[i], d, d_);
              io::imsave(d_, argv[i]);
            }
        }

      /// Compute the pw depth image from graph
      image2d<unsigned> maxdepth;
      resize(maxdepth, F);

      //const unsigned& (*maxptr) (const unsigned&, const unsigned&) = std::max<unsigned>;
      //write_vmap_to_image(g2, t, &tlink[0], gdepth, maxptr, 0u, maxdepth);

      image2d<value_t> rec;
      resize(rec, F);

      property_map<tree_t, V> vmap[NTREE];
      for (int i = 0; i < NTREE; ++i)
        vmap[i] = morpho::make_attribute_map_from_image(t[i],
                                                        imtransform(F, [i](value_t x) { return x[i]; }));
      auto colors = compute_graph_node_colors(g2, vmap);
      write_vmap_to_image_and_rec(g2, t, &tlink[0], gdepth, colors, std::greater<unsigned>(), maxdepth, rec);

      // Save
      io::imsave(imcast<uint16>(maxdepth), outname);
      io::imsave(rec, std::string(outname) + "-rec.tiff");
    }
  else if (attribute == OPT_COUNT)
    {
      image2d<unsigned> maxdepth = compute_map_count(g2, t,  &tlink[0]);

      // Save
      io::imsave(imcast<uint16>(maxdepth), outname);
    }
  else if (attribute == OPT_VARIATION)
    {
      /// Compute individual variation (for debugging)
      if (argc > 0)
        {
          image2d<unsigned> d_;
          resize(d_, F);

          for (int i = 0; i < std::min<int>(argc, NTREE); ++i)
            {
              auto vmap = morpho::make_attribute_map_from_image(t[i],
                                                                imtransform(F, [i](value_t x) { return x[i]; }));
              property_map<tree_t, uint16> var(t[i], 0);
              mln_foreach(auto x, t[i].nodes_without_root())
                var[x] = var[x.parent()] + std::abs(vmap[x] - vmap[x.parent()]);

              morpho::reconstruction(t[i], var, d_);
              io::imsave(d_, argv[i]);
            }
        }

      // Computer per-tree node value.
      property_map<tree_t, V> vmap[NTREE];
      for (int i = 0; i < NTREE; ++i)
        vmap[i] = morpho::make_attribute_map_from_image(t[i],
                                                        imtransform(F, [i](value_t x) { return x[i]; }));

      // Compute per-node colors
      auto colors = compute_graph_node_colors(g2, vmap);

      auto variation = compute_graph_variation(g2, colors);

      /// Compute the pw value image from graph
      image2d<unsigned> maxvar;
      resize(maxvar, F);

      const unsigned& (*maxptr) (const unsigned&, const unsigned&) = std::max<unsigned>;
      write_vmap_to_image(g2, t, &tlink[0], variation, maxptr, 0u, maxvar);

      // Save
      io::imsave(imcast<uint16>(maxvar), outname);
    }
  else if (attribute == OPT_AREA)
    {
      auto amap = compute_graph_area(g2, F, t, &tlink[0]);

      /// Compute the pw value image from graph
      image2d<float> maxvar;
      resize(maxvar, F);

      const float& (*maxptr) (const float&, const float&) = std::min<float>;
      write_vmap_to_image(g2, t, &tlink[0], amap, maxptr, 1.f, maxvar);

      // Save
      //io::imsave(maxvar, outname);
      io::imsave(imcast<uint16>(maxvar * (float)(1 << 16)), outname);
    }
}
