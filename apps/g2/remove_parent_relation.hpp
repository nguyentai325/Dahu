#ifndef REMOVE_PARENT_RELATION_HPP
# define REMOVE_PARENT_RELATION_HPP

# include <boost/property_map/vector_property_map.hpp>
# include <boost/graph/depth_first_search.hpp>
# include "types.hpp"

namespace mln
{

  /// G must be a tree !!!!
  template <class Graph>
  void
  update_parent_relation(Graph& g, tree_t& t1, tree_t& t2, tree_t& t3,
                         tlink_t& t1link, tlink_t& t2link, tlink_t& t3link)
  {

    typedef typename Graph::vertex_descriptor	vtype;
    typedef typename Graph::edge_descriptor	etype;

    typedef vec<vtype, 3> vec3v;
    boost::vector_property_map<vec3v> parent(boost::num_vertices(g));


    auto glink = boost::get(&graph_content::tlinks, g);

    //vtype root = boost::vertex(0, g);
    //parent[root] = vec3v{root, root, root};

    struct viz : public boost::default_dfs_visitor
    {
      viz(boost::vector_property_map<vec3v>& pmap)
        : parent(pmap)
      {
      }

      /// Initialize with the root node
      void initialize_vertex(vtype v, const Graph&)
      {
        parent[v] = vec3v{ v, v, v };
      }


      // void finish_vertex(vtype v, const Graph&)
      // {
      // 	(void) v;
      // 	std::cout << "Finish: " << v << " par: " << parent[v] << std::endl;
      // }


      void finish_edge(etype e, const Graph& g) const
      {
        auto glink = boost::get(&graph_content::tlinks, g);
        //auto gdepth = boost::get(&graph_content::depth, g);
        //auto gid = boost::get(boost::vertex_index, g);

        vtype s = boost::target(e, g); // parent
        vtype t = boost::source(e, g); // son

        for (int k = 0; k < 3; ++k)
          {
            if (glink[s][k].id() != tree_t::npos()) {
              parent[t][k] = s;
            } else {
              parent[t][k] = parent[s][k];
            }

            // else if (gid[parent[s][k]] == gid[s])
            //   parent[s][k] = parent[t][k];
            // else {
            //   unsigned d1 = gdepth[ parent[t][k] ][k];
            //   unsigned d2 = gdepth[ parent[s][k] ][k];
            //   if (d1 > d2)
            //     parent[s][k] = parent[t][k];
            // }
            //assert(gid[parent[t][k]] != gid[t]);
          }
        // std::cout << "Processing: " << e << std::endl;
      }


    private:
      boost::vector_property_map<vec3v>& parent;
    };

    // static_assert( boost::detail::has_member_function_finish_edge<viz, void>::value, "ya pas");
    // G is supposed to be a tree
    assert (boost::num_edges(g) == (boost::num_vertices(g)-1));

    std::cout << "Computing new parent" << std::endl;
    boost::depth_first_search(g, boost::visitor(viz(parent)) );
    std::cout << "End Computing new parent" << std::endl;

    // Update the senc of the graph
    {
      auto senc = boost::get(&graph_content::senc, g);
      BOOST_FOREACH(typename Graph::vertex_descriptor v, boost::vertices(g))
        {
          for (int k = 0; k < 3; ++k)
            senc[v][k] = (glink[v][k].id() != tree_t::npos()) ? glink[v][k].id() : glink[parent[v][k]][k].id();
        }
    }

    // Update now the trees
    tree_t* trees[3] = {&t1, &t2, &t3 };
    tlink_t* tlinks[3] = {&t1link, &t2link, &t3link };

    for (int k = 0; k < 3; ++k)
      {
        tlink_t& tlink = *(tlinks[k]);
        tree_t& tree = *(trees[k]);

        auto data = tree._get_data();

        // We need to process the node in reverse order
        // to update correctly the next_sibling tree
        mln_riter(cur, tree.nodes_without_root());
        cur.init();
        while (!cur.finished())
          {
            tree_t::node_type x = *cur;
            cur.next();

            vtype par = parent[tlink[x]][k];
            tree_t::node_type newpar = glink[par][k];

            // std::cout << tlink[x.id()] << " (" << x.id() << ")"
            //           << " par: " << tlink[x.parent().id()] << "(" << x.parent().id() << ")"
            //           << " newpar: " << par << " (" << newpar.id() << ")"
            //           << std::endl;
            //assert( k != 1 or newpar == x.parent() );
            if (newpar != x.parent())
              {
                tree_t::node_type lastnode = x.next_sibling().prev_node();

                // Remove the node (...)
                // 1. From the dble-linked list prev/next
                data->m_nodes[x.get_prev_node_id()].m_next = x.get_next_sibling_id();
                data->m_nodes[x.get_next_sibling_id()].m_prev = x.get_prev_node_id();
                // 2. From the next_sibling tree
                //std::cout << "Remove: Updating next sibling" << std::endl;
                for (tree_t::node_type y = x.prev_node(); y.next_sibling() == x; y = y.parent())
                  data->m_nodes[y.id()].m_next_sibling = x.get_next_sibling_id();


                // And insert the node at the right position
                // As a first child
                // 1. set parent
                data->m_nodes[x.id()].m_parent = newpar.id();
                // 2. insert in the dble-linked list prev/next
                data->m_nodes[lastnode.id()].m_next = newpar.get_next_node_id();
                data->m_nodes[newpar.get_next_node_id()].m_prev = lastnode.id();
                data->m_nodes[x.id()].m_prev = newpar.id();
                data->m_nodes[newpar.id()].m_next = x.id();
                // 3. update the next-sibling relation
                tree_t::vertex_id_t nexts = lastnode.get_next_node_id();
                //std::cout << "Insert: Updating next sibling" << std::endl;
                for (; lastnode != x; lastnode = lastnode.parent())
                  data->m_nodes[lastnode.id()].m_next_sibling = nexts;
              }
          }
      }
  }



}

#endif // ! REMOVE_PARENT_RELATION_HPP
