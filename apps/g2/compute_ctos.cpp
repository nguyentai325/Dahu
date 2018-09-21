#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/dontcare.hpp>
#include <mln/accu/accumulators/max.hpp>
#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/io/imsave.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include "compute_ctos.hpp"
#include "compute_g2.hpp"
#include "satmaxtree.hpp"

#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <tbb/parallel_for.h>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace mln
{

  template <class MyGraph>
  boost::vector_property_map<unsigned>
  compute_graph_depth(const MyGraph& g)
  {
    mln_entering("Compute graph depth");

    boost::vector_property_map<unsigned> depth(boost::num_vertices(g));

    auto one = [](mln::dontcare_t) -> int{ return 1; };
    auto w = boost::make_function_property_map<typename MyGraph::edge_descriptor, int, decltype(one)>(one);

    typename MyGraph::vertex_descriptor root = boost::vertex(0, g);
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
  compute_pw_count(const MyGraph& g,
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

  template <class MyGraph, class ValueMap, class BinaryFunction, class ValueType>
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
          typename MyGraph::vertex_descriptor gnode = tlink[i][tnode];
          w = op(w, vmap[gnode]);
        }
      px->val() = w;
    }

  }

  static
  void
  export_marginal_depth(const tree_t* t, int sz, std::string stem)
  {
    image2d<uint16> f;
    resize(f, t[0]._get_data()->m_pmap);

    for (int i = 0; i < sz; ++i)
      {
        auto depth = morpho::compute_depth(t[i]);
        morpho::reconstruction(t[i], depth, f);
        io::imsave(f, (boost::format("%s-%i.tiff") % stem % i).str());
      }
  }


  template <class V>
  tree_t
  compute_ctos(const image2d<V>& input,
               image2d<uint16>* imdepth,
               e_ctos_attribute attribute,
               ctos_extra_params_t params
               )
  {
    typedef V value_t;
    enum { NTREE = value_t::ndim };


    image2d<value_t> f = addborder_marginal(input);
    image2d<value_t> F = immerse_k1(f);
    image2d<uint16> maxdepth;
    point2d pmin = {0,0};
    resize(maxdepth, F);

    {
      /// Compute the marginal tree
      tree_t t[NTREE];
      tbb::parallel_for(0, (int)NTREE, [&t,&f,pmin](int i){
          t[i] = morpho::cToS_pinf(imtransform(f, [i](value_t x) {
                return x[i]; }), c4, pmin);
        });

      if (params.export_marginal_depth)
        export_marginal_depth(t, NTREE, params.export_marginal_depth_path);

      /// Compute the graph
      Graph<NTREE> g2;
      std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
      std::tie(g2, tlink) = compute_g2<NTREE>(t);

      if (attribute == CTOS_DEPTH)
        {
          /// Compute depth
          boost::vector_property_map<unsigned> gdepth;
          gdepth = compute_graph_depth(g2);


          /// Compute the pw depth image from graph
          write_vmap_to_image(g2, t, &tlink[0], gdepth,
                              functional::max_t<unsigned> (), uint16(0), maxdepth);
        }
      else if (attribute == CTOS_COUNT)
        {
          /// Compute depth
          boost::vector_property_map<unsigned> gdepth;
          gdepth = compute_graph_count(g2);


          /// Compute the pw depth image from graph
          write_vmap_to_image(g2, t, &tlink[0], gdepth,
                              functional::max_t<unsigned> (), uint16(0), maxdepth);

        }
      else if (attribute == CTOS_PW_COUNT)
        {
          auto maxdepth_ = compute_pw_count(g2, t,  &tlink[0]);
          unsigned maxv = accumulate(maxdepth_, accu::features::max<>());

          if (maxv > value_traits<uint16>::max())
            std::cerr << "Warning: too many values in the depth (depth image output is wrong)\n";

          maxdepth = eval(imcast<uint16>(maxdepth_));
        }
      else
        {
          std::cerr << "Unimplemented." << std::endl;
          std::abort();
        }
    }

    /// Compute the saturated maxtree
    tree_t tw;
    property_map<tree_t, uint16> depth2;
    std::tie(tw, depth2) = satmaxtree(maxdepth, pmin);

    /// Remove the 2F from the tree
    tw = tree_keep_2F(tw);

    /// Outputs the depth image
    if (imdepth) {
      resize(*imdepth, maxdepth);
      morpho::reconstruction(tw, depth2, *imdepth);
    }
    return tw;
  }


   tree_t compute_ctos_from_maxtrees(const image2d<rgb8>& input,
                                     image2d<uint16>* imdepth,
                                     bool mintree)
   {
     typedef rgb8 value_t;
     enum { NTREE = value_t::ndim };


     image2d<value_t> f = input; //addborder_marginal(input);
     image2d<value_t> F = input; //immerse_k1(f);


     /// Compute the marginal tree
     tree_t t[NTREE];
     tbb::parallel_for(0, (int)NTREE, [&t,&f,mintree](int i){
         if (mintree)
           t[i] = morpho::mintree_indexes(eval(channel(f, i)), c4);
         else
           t[i] = morpho::maxtree_indexes(eval(channel(f, i)), c4);
       });

     /// Compute the graph
     Graph<NTREE> g2;
     std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
     std::tie(g2, tlink) = compute_g2<NTREE>(t);


     /// Compute depth
     boost::vector_property_map<unsigned> gdepth;
     gdepth = compute_graph_depth(g2);


     /// Compute the pw depth image from graph
     image2d<uint16> maxdepth;
     resize(maxdepth, F);
     write_vmap_to_image(g2, t, &tlink[0], gdepth,
                         functional::max_t<unsigned> (), uint16(0), maxdepth);

     /// Compute the saturated maxtree
     tree_t tw;
     tw = morpho::maxtree_indexes(maxdepth, c4);

     /// Outputs the depth image
     if (imdepth) {
       *imdepth = maxdepth;
     }
     return tw;
   }

  // Explicit instanciation.
  template
  tree_t
  compute_ctos<rgb8>(const image2d<rgb8>& input,
                     image2d<uint16>* imdepth,
                     e_ctos_attribute attribute,
                     ctos_extra_params_t params
                     );

  // Explicit instanciation.
  template
  tree_t
  compute_ctos<rgb16>(const image2d<rgb16>& input,
                      image2d<uint16>* imdepth,
                      e_ctos_attribute attribute,
                      ctos_extra_params_t params
                      );

}
