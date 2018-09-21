#ifndef APPS_G2_COMPUTE_G2_HPP
# define APPS_G2_COMPUTE_G2_HPP

# include "types.hpp"

namespace mln
{
  /// \brief Compute the graph g2
  /// \returns the tuple (G, T₁, T₂, T₃) where
  ///          + G is the graph of shapes
  ///          + T₁ is the property map linking the first tree to a graph node
  ///          + T₂ is the property map linking the second tree to a graph node
  ///          + T₃ is the property map linking the third tree to a graph node
  // std::tuple<Graph,
  //            property_map<tree_t, Graph::vertex_descriptor>,
  //            property_map<tree_t, Graph::vertex_descriptor>,
  //            property_map<tree_t, Graph::vertex_descriptor> >
  // compute_g2(const tree_t& t1,
  //            const tree_t& t2,
  //            const tree_t& t3);


  /// \brief Generic version of the procedure to compute g2
  /// \param trees A vector of size NTREE containing the marginal trees
  template <unsigned NTREE>
  std::tuple<Graph<NTREE>, std::array<tlink_t, NTREE> >
  compute_g2(const tree_t* trees);


  // Explicit Declaration
  extern template
  std::tuple<Graph<2>, std::array<tlink_t, 2> >
  compute_g2<2>(const tree_t* trees);

  extern template
  std::tuple<Graph<3>, std::array<tlink_t, 3> >
  compute_g2<3>(const tree_t* trees);

  extern template
  std::tuple<Graph<4>, std::array<tlink_t, 4> >
  compute_g2<4>(const tree_t* trees);


  /***************************/
  /**  Implementation       **/
  /***************************/

  /// FWD
  property_map<tree_t, typename tree_t::node_type>
  smallest_enclosing_shape(const tree_t& t1,
			   const tree_t& t2,
			   const property_map<tree_t, unsigned>& d2);

  void
  compute_g2_precomputation(const tree_t* trees, int NTREE,
                            property_map<tree_t, unsigned>* d,
                            property_map<tree_t, tree_t::node_type>* SES);


  template <unsigned NTREE>
  std::tuple<Graph<NTREE>, std::array<tlink_t, NTREE> >
  compute_g2(const tree_t* trees)
  {
    mln_entering("mln::compute_g2");


    // 0. Compute the depth attribute on each tree and
    //    Compute the smallest enclosing shape each VS each tree
    // {
    typedef property_map<tree_t, unsigned> depth_attribute_t;
    depth_attribute_t d[NTREE];

    typedef property_map<tree_t, tree_t::node_type> tree_assoc_t;
    tree_assoc_t SES[NTREE][NTREE];

    compute_g2_precomputation(trees, NTREE, d, (tree_assoc_t*) SES);
    // }


    // 1. Compute the smallest enclosing shape each VS each tree
    // {
    // }

    // 2. Build the graph
    mln_entering("mln::compute_g2 - graph-construction - vertices");
    //for (unsigned i = 0; i < NTREE; ++i)
    //  std::cerr << "T" << i << ": " << trees[i].realsize() << std::endl;

    typedef Graph<NTREE> MyGraph;
    typedef graph_content<NTREE> my_graph_content;

    MyGraph graph(trees[0].realsize());

    // 2.1 Insert the nodes (not twice) and keep the correspondance tree <-> graph
    // For now, the graph is nor reduced, it does not matter.
    //{
    typedef property_map<tree_t, typename MyGraph::vertex_descriptor> graph_assoc_t;
    std::array<graph_assoc_t, NTREE> tlink;
    for (unsigned i = 0; i < NTREE; ++i)
      tlink[i] = graph_assoc_t (trees[i]);

    auto glink = boost::get(&my_graph_content::tlinks, graph);
    auto ulink = boost::get(&my_graph_content::ulink, graph);
    auto gdepth = boost::get(&my_graph_content::depth, graph);
    auto senc = boost::get(&my_graph_content::senc, graph);

    // Handle root seperatly
    typename MyGraph::vertex_descriptor v;
    typename MyGraph::vertex_iterator vcur, vend;

    // Root
    {
      std::tie(vcur, vend) = boost::vertices(graph);
      v = *vcur++;
      for (unsigned i = 0; i < NTREE; ++i)
        {
          tlink[i][trees[i].get_root_id()] = v;
          glink[v][i] = trees[i].get_root();
          gdepth[v][i] = 0;
          senc[v][i] = trees[i].get_root_id();
        }
      ulink[v] = trees[0].get_root();
    }

    // First tree : Special case since every node must be inserted
    mln_foreach(auto node, trees[0].nodes_without_root())
      {
        v = *vcur++;
        glink[v][0] = node; // Set graph -> tree link
        senc[v][0] = node.id();
        gdepth[v][0] = d[0][senc[v][0]];
        ulink[v] = node;

        for (unsigned i = 1; i < NTREE; ++i) {
          glink[v][i] = trees[i].nend();
          senc[v][i] = SES[0][i][node].id();
          gdepth[v][i] = d[i][senc[v][i]];
        }

        tlink[0][node] = v;   // Set tree -> graph link
      }

    // Other tree : Insert only if the node does not exist
    for (unsigned i = 1; i < NTREE; ++i)
      {
        const tree_t& t = trees[i];

        mln_foreach(auto node, t.nodes_without_root())
          {
            unsigned j = 0;
            for (j = 0; j < i; ++j)
              {
                tree_t::node_type ses_in_tj = SES[i][j][node];
                if (SES[j][i][ses_in_tj] == node) { // Already inserted
                  v = tlink[j][ses_in_tj];
                  break;
                }
              }
            if (j == i) // No node found, create new one
              {
                v = boost::add_vertex(graph);
                for (unsigned j = 0; j < NTREE; ++j)
                  if (i != j)
                    {
                      glink[v][j] = trees[j].nend();
                      senc[v][j] = SES[i][j][node].id();
                      gdepth[v][j] = d[j][senc[v][j]];
                    }
                ulink[v] = node;
              }
            senc[v][i] = node.id();
            gdepth[v][i] = d[i][node];
            glink[v][i] = node; // Set graph -> tree link
            tlink[i][node] = v;   // Set tree -> graph link
          }
      }
    mln_exiting();

    // Add the edges
    {
      mln_entering("mln::compute_g2 - graph-construction - edges");
      typename MyGraph::vertex_iterator cur, end;
      boost::tie(cur,end) = boost::vertices(graph);
      for (++cur; cur != end; ++cur) //skip the root
        {
          for (unsigned i = 0; i < NTREE; ++i) {
            tree_t::vertex_id_t par;
            tree_t::node_type s = glink[*cur][i];
            if (s != trees[i].nend())
              par = s.get_parent_id();
            else
              par = senc[*cur][i];
            boost::add_edge(*cur, tlink[i][par], graph);
          }
        }
      mln_exiting();
    }
    //}

    // 3. Reduction step
    {
      typename MyGraph::vertex_iterator v, vend;
      typename MyGraph::out_edge_iterator e1, e2, next, eend;
      boost::tie(v,vend) = boost::vertices(graph);
      for (; v != vend; ++v)
      	{
      	  std::tie(e1, eend) = boost::out_edges(*v, graph);
      	  for (next = e1; e1 != eend; e1 = next)
      	    {
      	      typename MyGraph::vertex_descriptor v1 = boost::target(*e1, graph);
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
    //std::cerr << "Graph: " << boost::num_vertices(graph) << std::endl;

    mln_exiting();
    return std::tie(graph, tlink);
  }

}

#endif // ! APPS_G2_COMPUTE_G2_HPP
