#ifndef APPS_G2_ROUTINES_HPP
# define APPS_G2_ROUTINES_HPP

# include "types.hpp"
# include <boost/property_map/vector_property_map.hpp>

# include <mln/core/colors.hpp>
# include <mln/core/grays.hpp>

namespace mln
{

  ///
  /// \brief Compute a per-node color for g2
  ///
  /// The color for a node \p x is `[rX, gY, bZ]` where
  /// + rX is the value of the smallest enclosing red shape X of x
  /// + gY is the value of the smallest enclosing green shape Y of x
  /// + bZ is the value of the smallest enclosing blue shape Z of x
  template <typename V>
  boost::vector_property_map< vec<V, NTREE> >
  compute_graph_node_colors(const MyGraph& g2,
                            const property_map<tree_t, V>* vmap);

                            // const property_map<tree_t, uint8>& vr,
                            // const property_map<tree_t, uint8>& vg,
                            // const property_map<tree_t, uint8>& vb);

  ///
  /// \brief compute the indicatrice for each point of the LB of the nodes
  /// in the graph.
  ///
  /// I(x) = <Rx, Gx, Bx>
  /// where Rx, Gx, Bx are graph node vertex ids or (-1) if the node ∉ LB
  ///
  ///
  /// with:
  ///  S¹ = SES(red, {x})
  ///  S² = SES(green, {x})
  ///  S³ = SES(blue, {x})
  ///  LB = LOWER_BOUNDS_⊆({S¹, S², S³})
  ///  Rx = S¹ if S¹ ∈ LB else ∅
  ///  Gx = S² if S² ∈ LB else ∅
  ///  Bx = S³ if S³ ∈ LB else ∅
  void
  compute_graph_map(const MyGraph& g2,
                    const tree_t& t1, const tree_t& t2, const tree_t& t3,
                    const tlink_t& t1link, const tlink_t& t2link, const tlink_t& t3link,
                    image2d<vec3u>& out);


  void
  compute_graph_map(const MyGraph& g2,
                    const tree_t* trees,
                    const tlink_t* tlink,
                    image2d< vec<unsigned, NTREE> >& out);



  /************************************/
  /** Implementation               ****/
  /************************************/

  template <typename V>
  boost::vector_property_map< vec<V,NTREE> >
  compute_graph_node_colors(const MyGraph& g2,
                            const property_map<tree_t, V>* vmap)
  {
    mln_entering("mln::compute_graph_node_colors");


    boost::vector_property_map< vec<V,NTREE> > output(boost::num_vertices(g2));
    auto senc = boost::get(&my_graph_content::senc, g2);

    MyGraph::vertex_iterator vi, vend;
    boost::tie(vi, vend) = boost::vertices(g2);
    for (; vi != vend; ++vi)
      for (int k = 0; k < NTREE; ++k)
        {
          unsigned tnode = senc[*vi][k];
          output[*vi][k] = vmap[k][tnode];
        }

    mln_exiting();
    return output;
  }

}

#endif // ! APPS_G2_ROUTINES_HPP
