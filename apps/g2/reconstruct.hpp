#ifndef RECONSTRUCT_HPP
# define RECONSTRUCT_HPP

# include <mln/morpho/component_tree/filtering.hpp>
# include <boost/foreach.hpp>
# include "types.hpp"

namespace mln
{

  template <class Graph, class AMap>
  void
  reconstruct_marginal(const tree_t& t1, const tree_t& t2, const tree_t& t3,
		       const AMap& val1, const AMap& val2, const AMap& val3,
		       const Graph& graph,
		       image2d<rgb8>& out);

  template <class Graph, class ColorMap>
  image2d<rgb8>
  reconstruct_from_graph(const Graph& g, const image2d<vec3u>& pmap,
                         const ColorMap& vmap);



  /***************************/
  /**   Implementation      **/
  /***************************/


  template <class Graph, class AMap>
  void
  reconstruct_marginal(const tree_t& t1, const tree_t& t2, const tree_t& t3,
		       const AMap& val1, const AMap& val2, const AMap& val3,
		       const Graph& graph,
		       image2d<rgb8>& out)
  {
    auto glinks = boost::get(&graph_content::tlinks, graph);
    const tree_t* t_array[3] = {&t1, &t2, &t3};
    //const AMap* val_array[3] = {&val1, &val2, &val3};

    // 1. Get the active nodes in the trees from graph nodes.
    property_map<tree_t, bool> active[3] = {
      property_map<tree_t, bool> (t1, false),
      property_map<tree_t, bool> (t2, false),
      property_map<tree_t, bool> (t3, false)
    };

    BOOST_FOREACH(typename Graph::vertex_descriptor v, boost::vertices(graph))
      {
        for (int i = 0; i < 3; ++i)
          if (glinks[v][i] != t_array[i]->nend() )
            active[i][glinks[v][i]] = true;
      }

    filter_direct_and_reconstruct(t1, active[0], val1, red(out));
    filter_direct_and_reconstruct(t2, active[1], val2, green(out));
    filter_direct_and_reconstruct(t3, active[2], val3, blue(out));
  }

  template <class graph_t, class AMap>
  image2d<rgb8>
  reconstruct_from_graph(const Graph& g, const graph_t& fg,
                         const tree_t& t1, const tree_t& t2, const tree_t& t3,
                         const AMap& val1, const AMap& val2, const AMap& val3,
                         const image2d<vec3u>& pmap)

  {
    // 1. Reconstruct the per-node new values.
    auto glinks = boost::get(&graph_content::tlinks, g);
    const tree_t* t_array[3] = {&t1, &t2, &t3};


    property_map<tree_t, bool> active[3] = {
      property_map<tree_t, bool> (t1, false),
      property_map<tree_t, bool> (t2, false),
      property_map<tree_t, bool> (t3, false)
    };

    BOOST_FOREACH(typename Graph::vertex_descriptor v, boost::vertices(fg))
      {
        for (int i = 0; i < 3; ++i)
          if (glinks[v][i] != t_array[i]->nend() )
            active[i][glinks[v][i]] = true;
      }

    auto vmap1 = filter_direct(t1, active[0], val1);
    auto vmap2 = filter_direct(t2, active[1], val2);
    auto vmap3 = filter_direct(t3, active[2], val3);


    auto vmap = compute_graph_node_colors(g, vmap1, vmap2, vmap3);
    // BOOST_FOREACH(typename Graph::vertex_descriptor v, boost::vertices(g))
    //   std::cout << "node " << v << " : " << (vec3u)vmap[v] << std::endl;
    // std::cout << "=========" << std::endl;

    // 2. reconstruct from the graph
    image2d<rgb8> out;
    resize(out, pmap);

    mln_pixter(pxin, pxout, pmap, out);
    mln_forall(pxin, pxout)
      {
        vec3u p = pxin->val();
        rgb8 v = { vmap[p[0]][0],
                   vmap[p[1]][1],
                   vmap[p[2]][2] };

        // if (p[0] != (unsigned) -1) {
        //   v[0] = vmap[p[0]][0];
        //   if (p[1] == (unsigned)-1) v[1] = vmap[p[0]][1];
        //   if (p[2] == (unsigned)-1) v[2] = vmap[p[0]][2];
        // }
        // if (p[1] != (unsigned) -1) {
        //   v[1] = vmap[p[1]][1];
        //   if (p[0] == (unsigned)-1) v[0] = vmap[p[1]][0];
        //   if (p[2] == (unsigned)-1) v[2] = vmap[p[1]][2];
        // }
        // if (p[2] != (unsigned) -1) {
        //   v[2] = vmap[p[2]][2];
        //   if (p[0] == (unsigned)-1) v[0] = vmap[p[2]][0];
        //   if (p[1] == (unsigned)-1) v[1] = vmap[p[2]][1];
        // }

        pxout->val() = v;
      }
    return out;
  }


}

#endif // ! RECONSTRUCT_HPP
