#include "routines.hpp"
#include <mln/core/trace.hpp>

namespace mln
{

  void
  compute_graph_map(const MyGraph& g2,
                    const tree_t& t1, const tree_t& t2, const tree_t& t3,
                    const tlink_t& t1link, const tlink_t& t2link, const tlink_t& t3link,
                    image2d<vec3u>& out)
  {
    mln_entering("mln::compute_graph_map");

    auto SES = boost::get(&my_graph_content::senc, g2);

    mln_pixter(px, out);
    mln_forall(px)
      {
        tree_t::vertex_id_t n[3];
        n[0] = t1.get_node_id(px->index());
        n[1] = t2.get_node_id(px->index());
        n[2] = t3.get_node_id(px->index());

        // Compute the LB of {n1, n2, n3}
        unsigned gnode[3];
        gnode[0] = t1link[n[0]];
        gnode[1] = t2link[n[1]];
        gnode[2] = t3link[n[2]];

        auto inc = [n, &gnode, &SES](int a, int b) {
          return SES[gnode[a]][b] == n[b];
        };


        if (inc(1,0)) gnode[0] = gnode[1];
        if (inc(2,0)) gnode[0] = gnode[2];
        if (inc(0,1)) gnode[1] = gnode[0];
        if (inc(2,1)) gnode[1] = gnode[2];
        if (inc(0,2)) gnode[2] = gnode[0];
        if (inc(1,2)) gnode[2] = gnode[1];

        px->val() = {gnode[0], gnode[1], gnode[2]};
      }

    mln_exiting();
  }


  void
  compute_graph_map(const MyGraph& g2,
                    const tree_t* t,
                    const tlink_t* tlink,
                    image2d< vec<unsigned, NTREE> >& out)
  {
    mln_entering("mln::compute_graph_map");

    auto d = boost::get(&my_graph_content::depth, g2);

    productorder_less< vec<unsigned, NTREE> > prec;

    mln_pixter(px, out);
    mln_forall(px)
      {
        unsigned gnode[NTREE];

        for (int i = 0; i < NTREE; ++i)
          {
            tree_t::vertex_id_t n = t[i].get_node_id(px->index());
            gnode[i] = tlink[i][n];
          }

        for (int i = 0; i < NTREE; ++i)
          {
            for (int j = 0; j < NTREE; ++j)
              if (i != j and prec(d[gnode[j]], d[gnode[i]]))
                gnode[i] = gnode[j];
            px->val()[i] = gnode[i];
          }
      }

    mln_exiting();
  }

}


