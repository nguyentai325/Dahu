#include "compute_g2.hpp"

# include <mln/morpho/component_tree/accumulate.hpp>
# include <mln/morpho/component_tree/compute_depth.hpp>
# include <apps/g2/accu/lca.hpp>
# include <tbb/parallel_for.h>

namespace mln
{

  property_map<tree_t, typename tree_t::node_type>
  smallest_enclosing_shape(const tree_t& t1,
			   const tree_t& t2,
			   const property_map<tree_t, unsigned>& d2)
  {
    mln_entering("smallest_enclosing_shape");
    accu::least_common_ancestor<unsigned, image2d<unsigned> >  acc(t2, d2);
    auto tmp = morpho::accumulate(t1, acc);
    mln_exiting();
    return tmp;
  }


  void
  compute_g2_precomputation(const tree_t* trees, int NTREE,
                            property_map<tree_t, unsigned>* d,
                            property_map<tree_t, tree_t::node_type>* SES)
  {
    mln_entering("G2 - precomputation");

    // 0. Compute the depth attribute on each tree
    // {
    for (int i = 0; i < NTREE; ++i)
      d[i] = morpho::compute_depth(trees[i]);

    // }

    // 1. Compute the smallest enclosing shape each VS each tree
    // {
    tbb::parallel_for(0, (int)NTREE, [&SES,&trees,&d,NTREE](int i) {
        for (int j = 0; j < NTREE; ++j)
          if (i != j)
            SES[i * NTREE + j] = smallest_enclosing_shape(trees[i], trees[j], d[j]);
      });
    // }

    mln_exiting();
  }



  /*
  std::tuple<Graph,
	     property_map<tree_t, Graph::vertex_descriptor>,
	     property_map<tree_t, Graph::vertex_descriptor>,
	     property_map<tree_t, Graph::vertex_descriptor> >
  compute_g2(const tree_t& t1,
             const tree_t& t2,
             const tree_t& t3)
  {
    mln_entering("mln::compute_g2");

    // typedef uint8 V;

    // image2d<V> r = transform(f, [](rgb8 x) -> V { return x[0]; });
    // image2d<V> g = transform(f, [](rgb8 x) -> V { return x[1]; });
    // image2d<V> b = transform(f, [](rgb8 x) -> V { return x[2]; });

    auto d1 = morpho::compute_depth(t1);
    auto d2 = morpho::compute_depth(t2);
    auto d3 = morpho::compute_depth(t3);


    // 1. Compute the smallest enclosing shape
    // {
    auto SES_12 = smallest_enclosing_shape(t1, t2, d2);
    auto SES_13 = smallest_enclosing_shape(t1, t3, d3);
    auto SES_21 = smallest_enclosing_shape(t2, t1, d1);
    auto SES_23 = smallest_enclosing_shape(t2, t3, d3);
    auto SES_31 = smallest_enclosing_shape(t3, t1, d1);
    auto SES_32 = smallest_enclosing_shape(t3, t2, d2);
    // }

    // 2. Build the graph
    mln_entering("mln::compute_g2 - graph-construction - vertices");

    Graph graph(t1.size()-1);

    // 2.1 Insert the nodes (not twice) and keep the correspondance tree <-> graph
    // For now, the graph is nor reduced, it does not matter.
    //{
    property_map<tree_t, Graph::vertex_descriptor> t1link(t1);
    property_map<tree_t, Graph::vertex_descriptor> t2link(t2);
    property_map<tree_t, Graph::vertex_descriptor> t3link(t3);

    auto glink = boost::get(&graph_content::tlinks, graph);
    auto ulink = boost::get(&graph_content::ulink, graph);
    auto gdepth = boost::get(&graph_content::depth, graph);
    auto senc = boost::get(&graph_content::senc, graph);

    // Handle root seperatly
    Graph::vertex_descriptor v;
    auto vcurrent = boost::vertices(graph).first;
    v = *vcurrent++;

    t1link[t1.get_root_id()] = v;
    t2link[t2.get_root_id()] = v;
    t3link[t3.get_root_id()] = v;
    glink[v][0] = t1.get_root();
    glink[v][1] = t2.get_root();
    glink[v][2] = t3.get_root();
    gdepth[v][0] = 0;
    gdepth[v][1] = 0;
    gdepth[v][2] = 0;
    senc[v][0] = t1.get_root_id();
    senc[v][1] = t2.get_root_id();
    senc[v][2] = t3.get_root_id();
    ulink[v] = t1.get_root();

    // First tree
    mln_foreach(auto node, t1.nodes_without_root())
      {
        v = *vcurrent;
        glink[v][0] = node; // Set graph -> tree link
        glink[v][1] = t2.nend();
        glink[v][2] = t3.nend();
        ulink[v] = node;
        senc[v][0] = node.id();
        senc[v][1] = SES_12[node].id();
        senc[v][2] = SES_13[node].id();
        gdepth[v][0] = d1[senc[v][0]];
        gdepth[v][1] = d2[senc[v][1]];
        gdepth[v][2] = d3[senc[v][2]];

        t1link[node] = v;   // Set tree -> graph link
        ++vcurrent;
      }

    // Second tree
    mln_foreach(auto node, t2.nodes_without_root())
      {
        if (SES_12[SES_21[node]] == node) // Already inserted
          v = t1link[SES_21[node]];
        else {
          v = boost::add_vertex(graph);
          glink[v][0] = t1.nend();
          glink[v][2] = t3.nend();
          senc[v][0] = SES_21[node].id();
          senc[v][1] = node.id();
          senc[v][2] = SES_23[node].id();
          gdepth[v][0] = d1[senc[v][0]];
          gdepth[v][1] = d2[senc[v][1]];
          gdepth[v][2] = d3[senc[v][2]];
          ulink[v] = node;
        }
        glink[v][1] = node; // Set graph -> tree link
        t2link[node] = v;   // Set tree -> graph link
      }

    // Third tree
    mln_foreach(auto node, t3.nodes_without_root())
      {
        if (SES_13[SES_31[node]] == node) // Already inserted
          v = t1link[SES_31[node]];
        else if (SES_23[SES_32[node]] == node) // Already inserted
          v = t2link[SES_32[node]];
        else {
          v = boost::add_vertex(graph);
          glink[v][0] = t1.nend();
          glink[v][1] = t2.nend();
          senc[v][0] = SES_31[node].id();
          senc[v][1] = SES_32[node].id();
          senc[v][2] = node.id();
          gdepth[v][0] = d1[senc[v][0]];
          gdepth[v][1] = d2[senc[v][1]];
          gdepth[v][2] = d3[senc[v][2]];
          ulink[v] = node;
        }
        glink[v][2] = node; // Set graph -> tree link
        t3link[node] = v;   // Set tree -> graph link
      }
    mln_exiting();

    // Add the edges
    {
      mln_entering("mln::compute_g2 - graph-construction - edges")
      Graph::vertex_iterator cur, end;
      boost::tie(cur,end) = boost::vertices(graph);
      for (++cur; cur != end; ++cur) //skip the root
        {
          tree_t::node_type s1 = glink[*cur][0];
          tree_t::node_type s2 = glink[*cur][1];
          tree_t::node_type s3 = glink[*cur][2];

          bool has_s1 = (s1 != t1.nend());
          bool has_s2 = (s2 != t2.nend());
          bool has_s3 = (s3 != t3.nend());

          tree_t::node_type p1 = has_s1 ? s1.parent() : (has_s2 ? SES_21[s2] : SES_31[s3]);
          tree_t::node_type p2 = has_s2 ? s2.parent() : (has_s1 ? SES_12[s1] : SES_32[s3]);
          tree_t::node_type p3 = has_s3 ? s3.parent() : (has_s1 ? SES_13[s1] : SES_23[s2]);

          boost::add_edge(*cur, t1link[p1], graph);
          boost::add_edge(*cur, t2link[p2], graph);
          boost::add_edge(*cur, t3link[p3], graph);
        }
      mln_exiting();
    }
    //}

    // 3. Reduction step
    {
      Graph::vertex_iterator v, vend;
      Graph::out_edge_iterator e1, e2, next, eend;
      boost::tie(v,vend) = boost::vertices(graph);
      for (; v != vend; ++v)
      	{
      	  std::tie(e1, eend) = boost::out_edges(*v, graph);
      	  for (next = e1; e1 != eend; e1 = next)
      	    {
      	      Graph::vertex_descriptor v1 = boost::target(*e1, graph);
      	      ++next;
      	      e2 = boost::out_edges(*v, graph).first;
      	      for (; e2 != eend; ++e2)
                if (vecprod_isless(gdepth[v1], gdepth[boost::target(*e2, graph)]))
                  {
                    //std::cout << "Remove: " << *e1 << std::endl;
                    boost::remove_edge(*e1, graph);
                    break;
                  }
      	    }
      	  // std::tie(e1, eend) = boost::out_edges(*v, graph);
	  // std::cout << "--";
	  // for (; e1 != eend; ++e1)
	  //    std::cout << *e1;
	  // std::cout << std::endl;
      	}
    }



    // Some stat
    std::cout << "T1: " << t1.size() << std::endl;
    std::cout << "T2: " << t2.size() << std::endl;
    std::cout << "T3: " << t3.size() << std::endl;
    std::cout << "Graph: " << boost::num_vertices(graph) << std::endl;

    mln_exiting();
    return std::tie(graph, t1link, t2link, t3link);
  }
  */


 // Explicit definition
 template
 std::tuple<Graph<2>, std::array<tlink_t, 2> >
 compute_g2<2>(const tree_t* trees);

 template
 std::tuple<Graph<3>, std::array<tlink_t, 3> >
 compute_g2<3>(const tree_t* trees);

 template
 std::tuple<Graph<4>, std::array<tlink_t, 4> >
 compute_g2<4>(const tree_t* trees);

}
