#ifndef MLN_MORPHO_COMPONENT_TREE_CUTS_HPP
# define MLN_MORPHO_COMPONENT_TREE_CUTS_HPP

/// \file
/// \brief This file provides some routines to process hirarchical segmentation
/// such as cuts on hierarchy

namespace mln
{

  namespace morpho
  {

    /// \brief Cut on hierarchy
    ///
    /// Given a decreasing criterion 𝓒, a component Γ is kept in the
    /// hierarchy if 𝓒(Γ) but ¬𝓒(parent(Γ))
    template <class P, class Amap, class CriterionMap, class ValueMap, class OutputImage>
    void
    cut_and_reconstruct(const component_tree<P, Amap>& tree,
                        const CriterionMap& criterion_map,
                        const ValueMap& value_map,
                        OutputImage&& out);

    template <class P, class Amap, class CriterionMap>
    void
    cut_inplace(component_tree<P, Amap>& tree,
                const CriterionMap& criterion_map);




    /******************************/
    /***   Implementation      ****/
    /******************************/

    template <class P, class Amap, class CriterionMap, class ValueMap, class OutputImage>
    void
    cut_and_reconstruct(const component_tree<P, Amap>& tree,
                        const CriterionMap& pred,
                        const ValueMap& value_map,
                        OutputImage&& out)
    {
      mln_entering("mln::morpho::cut_and_reconstruct");

      typedef component_tree<P, Amap> tree_t;
      typedef typename tree_t::vertex_id_t vertex_id_t;
      typedef typename ValueMap::value_type V;


      property_map<tree_t, V>           vmap(tree);
      vertex_id_t r = tree.get_root_id();
      vmap[r] = value_map[r];

      mln_foreach(auto x, tree.nodes_without_root())
        if (pred[x.parent()])
          {
            vmap[x] = vmap[x.parent()];
          }
        else
          {
            vmap[x] = value_map[x];
          }

      mln_foreach(auto px, out.pixels())
        {
          auto x = tree.get_node_at(px.index());
          px.val() = vmap[x];
        }

      mln_exiting();
    }

    template <class P, class Amap, class CriterionMap>
    void
    cut_inplace(component_tree<P, Amap>& tree,
                const CriterionMap& pred)
    {
      typedef component_tree<P, Amap> tree_t;
      typedef typename tree_t::vertex_id_t vertex_id_t;

      mln_entering("mln::morpho::cut_inplace");
      auto data = tree._get_data();

      // forward update parent
      for (auto x = tree.get_root(); x.id() != tree.npos();)
        {
          if (pred[x.parent()])
            {
              vertex_id_t q = x.get_parent_id();

              // unactive subtree and set parent
              auto end = x.next_sibling();
              for (; x != end; x = x.next_node()) {
                data->m_nodes[x.id()].m_parent = q;
              }
            }
          else
            {
              x = x.next_node();
            }
        }

      // backward update next-sibling and the dble-linked list
      mln_reverse_foreach(auto x, tree.nodes())
        {
          typename tree_t::node_type snext = x.next_sibling();
          if (pred[snext.parent()]) { // snext is removed
            vertex_id_t snextid = data->m_nodes[snext.id()].m_next;
            data->m_nodes[x.id()].m_next_sibling = snextid;
          }

          if (pred[x.parent()]) { // the current node is removed
            data->m_nodes[x.get_prev_node_id()].m_next = x.get_next_node_id();
            data->m_nodes[x.get_next_node_id()].m_prev = x.get_prev_node_id();
          }

          assert(not pred[x.parent().parent()]);
          assert(not pred[x.next_node().parent()]);
          assert(not pred[x.next_sibling().parent()]);
        }

      // reaffect point to nodes
      mln_foreach (auto& v, data->m_pmap.values())
        {
          if (pred[tree.get_node(v).parent()])
            v = data->m_nodes[v].m_parent;
        }
      mln_exiting();
    }


  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_CUTS_HPP
