#include <mln/core/neighb2d.hpp>
#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include "satmaxtree.hpp"

namespace mln
{

  std::pair<
    morpho::component_tree<unsigned, image2d<unsigned> >,
    property_map<morpho::component_tree<unsigned, image2d<unsigned> >, uint16>
    >
  satmaxtree(const image2d<uint16>& f, point2d pmin)
  {
    mln_entering("mln::satmaxtree");

    typedef morpho::component_tree<unsigned, image2d<unsigned> > T;

    // bool ASK = true;
    // point2d pmin{0,0};
    // if (ASK) {
    //   std::cout << "Point a l'infini: " << std::endl;
    //   std::cin >> pmin[0] >> pmin[1];
    //   pmin *= 2;
    // }
    T tree = morpho::cToS_pinf(f, c4, pmin);

    auto newdata1 = tree._get_data();



    box2d olddom1 = newdata1->m_pmap.domain();
    std::cout << "s11  "  << olddom1 << std::endl;

	std::cout << "minh dai ca vl"  << std::endl;
    property_map<T, uint16> vmap;
    {
      image2d<uint16> F = immerse_k1(f, 69);
      vmap = morpho::make_attribute_map_from_image(tree, F);
      vmap[tree.npos()] = 0;
          std::cout << "domain  "<< F.domain() << std::endl;

    }

    auto predfun = [&vmap,&tree](const T::vertex_id_t& n) {
      return vmap[n] > vmap[tree.get_node(n).parent()] or tree.get_node(n).get_parent_id() == tree.npos();
    };
    auto pred = make_functional_property_map<T::vertex_id_t> (predfun);

    // Subtractive
    {
      property_map<T, unsigned> delta(tree, 0);
      unsigned n = 0;
      mln_foreach(auto x, tree.nodes_without_root()) {
        delta[x] = delta[x.parent()];
        if (vmap[x] < vmap[x.parent()])
          delta[x] += vmap[x.parent()] - vmap[x];
        ++n;
      }
      // We need to separate the two loops !
      mln_foreach(auto x, tree.nodes_without_root()) {
        vmap[x] += delta[x];
      }
      //std::cerr << "Number of nodes before: " << n << std::endl;
    }

    morpho::filter_direct_inplace(tree, pred);
    {
      unsigned n = 0;
      mln_foreach(auto x, tree.nodes()) {
        (void) x;
        assert(vmap[x] > vmap[x.parent()]);
        ++n;
      }
      //std::cerr << "Number of nodes after: " << n << std::endl;
    }

    mln_exiting();
    return {std::move(tree), std::move(vmap)};
  }




  morpho::component_tree<unsigned, image2d<unsigned> >
  tree_keep_2F(const morpho::component_tree<unsigned, image2d<unsigned> >& tree)
  {
    typedef unsigned P;
    morpho::component_tree<P, image2d<P> > out;

    auto newdata = out._get_data();
    auto olddata = tree._get_data();

    // 1. Copy the point2node map
    box2d olddom = olddata->m_pmap.domain();
    box2d dom;
    dom.pmin = olddom.pmin / 2;
    dom.pmax = (olddom.pmax + 1) / 2;
    newdata->m_pmap.resize(dom);
    copy(olddata->m_pmap | sbox2d{olddom.pmin, olddom.pmax, {2,2}},
         newdata->m_pmap);

    // 2. Copy the node
    newdata->m_nodes = olddata->m_nodes;

    // 3. Copy the point set and update node first point index/
    newdata->m_S.resize(dom.size());
    unsigned j = 0;
    for (unsigned i = 0; i < olddata->m_S.size(); ++i)
      {
        P p = olddata->m_S[i];
        point2d pt = olddata->m_pmap.point_at_index(p);
        if (K1::is_face_2(pt))
          {
            newdata->m_S[j] = newdata->m_pmap.index_of_point(pt/2);
            auto node = tree.get_node_at(p);
            if (node.get_first_point_id() == i)
              newdata->m_nodes[node.id()].m_point_index = j;
            ++j;
          }
      }
    // 4. Do not forget the sentinel
    newdata->m_nodes[out.npos()].m_point_index = j;

    return out.get_subtree(tree.get_root_id());
  }


}
