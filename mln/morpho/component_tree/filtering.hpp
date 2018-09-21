#ifndef MLN_MORPHO_COMPONENT_TREE_FILTERING_HPP
# define MLN_MORPHO_COMPONENT_TREE_FILTERING_HPP

# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace morpho
  {

    /// \brief Min filtering strategy
    ///
    /// A node is remove if it does not pass the predicate and
    /// all its descendant are removed as well.
    /// Formally, a component Γ is kept if 𝓒(Γ) with
    /// 𝓒(Γ) = ⋁ { 𝓟(Γ'), Γ ⊆ Γ' }.
    template <class P, class Amap, class PredicateMap>
    void
    filter_min_inplace(component_tree<P, Amap>& tree,
                       const PredicateMap& predicate_map);

    /// \brief Max filtering strategy
    ///
    /// A node is remove if every node in its childhood does not pass
    /// the predicate. All its descendant are removed as well.  Formally, a
    /// component Γ is kept if 𝓒(Γ) with 𝓒(Γ) = ⋀ { 𝓟(Γ'), Γ' ⊆ Γ }.
    template <class P, class Amap, class PredicateMap>
    void
    filter_max_inplace(component_tree<P, Amap>& tree,
                       const PredicateMap& predicate_map);


    /// \brief Direct filtering strategy
    /// A node is remove if it does not pass the predicate, the others
    /// remain.
    template <class P, class Amap, class PredicateMap>
    void
    filter_direct_inplace(component_tree<P, Amap>& tree,
                          const PredicateMap& predicate_map);


    /// \brief Direct filtering strategy
    /// A node is remove if it does not pass the predicate, the others
    /// remain.
    template <class P, class Amap, class PredicateMap, class ValueMap>
    property_map<component_tree<P, Amap>, typename ValueMap::value_type>
    filter_direct(const component_tree<P, Amap>& tree,
                  const PredicateMap& predicate_map,
                  const ValueMap& value_map);

    template <class P, class Amap, class PredicateMap, class ValueMap, class OutputImage>
    void
    filter_direct_and_reconstruct(const component_tree<P, Amap>& tree,
                                  const PredicateMap& predicate_map,
                                  const ValueMap& value_map,
                                  OutputImage&& out);



    /********************/
    /** Implementation **/
    /********************/

    namespace internal
    {

      /// \brief Prune a tree inplace.
      /// The criterion is supposed increasing.
      template <class P, class Amap, class Criterion>
      void
      prune_inplace(component_tree<P, Amap>& tree,
                    const Criterion& pred)
      {
        typedef component_tree<P, Amap> tree_t;
        typedef typename tree_t::vertex_id_t vertex_id_t;

        auto data = tree._get_data();

        property_map<tree_t, bool> alive(tree, true);

        // forward: update parent
        for (auto x = tree.get_root(); x.id() != tree.npos();)
          {
            if (not pred[x])
              {
                vertex_id_t q = x.get_parent_id();

                // unactive subtree and set parent
                auto end = x.next_sibling();
                for (; x != end; x = x.next_node()) {
                  alive[x] = false;
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
            vertex_id_t snext = x.get_next_sibling_id();
            if (not alive[snext]) {
              snext = data->m_nodes[snext].m_next;
              data->m_nodes[x.id()].m_next_sibling = snext;
            }

            if (not alive[x.id()]) {
              data->m_nodes[x.get_prev_node_id()].m_next = x.get_next_node_id();
              data->m_nodes[x.get_next_node_id()].m_prev = x.get_prev_node_id();
            }

            assert(alive[x.get_parent_id()]);
            assert(alive[x.get_next_node_id()]);
            assert(alive[x.get_next_sibling_id()]);
          }

        // reaffect point to nodes
        mln_foreach (auto& v, data->m_pmap.values())
          {
            if (not alive[v])
              v = data->m_nodes[v].m_parent;
          }
      }
    }


    template <class P, class Amap, class PredicateMap>
    void
    filter_min_inplace(component_tree<P, Amap>& tree,
                       const PredicateMap& pred)
    {
      mln_entering("mln::morpho::filter_min_inplace");

      // even if pred is not increasing, as
      // soon as pred[n] is false, the subtree is removed.
      internal::prune_inplace(tree, pred);
      mln_exiting();
    }

    template <class P, class Amap, class PredicateMap>
    void
    filter_max_inplace(component_tree<P, Amap>& tree,
                       const PredicateMap& pred)
    {
      mln_entering("mln::morpho::filter_max_inplace");

      typedef component_tree<P, Amap> tree_t;
      auto data = tree._get_data();

      property_map<tree_t, bool> active(tree, true);
      mln_reverse_foreach(const typename tree_t::node_type& n, tree.nodes())
        {
          active[n] = active[n] or pred[n];
          active[n.parent()] = active[n.parent()] or active[n];
        }

      internal::prune_inplace(tree, active);
      mln_exiting();
    }



    template <class P, class Amap, class PredicateMap>
    void
    filter_direct_inplace(component_tree<P, Amap>& tree,
                          const PredicateMap& pred)
    {
      mln_entering("mln::morpho::filter_direct_inplace");

      typedef component_tree<P, Amap> tree_t;
      typedef typename tree_t::vertex_id_t vertex_id_t;

      auto data = tree._get_data();

      // Forward: set parent
      for (auto x = tree.get_root(); x.id() != tree.npos(); x = x.next_node())
        {
          vertex_id_t q = x.get_parent_id();
          if (not pred[q])
            data->m_nodes[x.id()].m_parent = data->m_nodes[q].m_parent;
        }

      // Backward: set other relation
      mln_reverse_foreach(auto x, tree.nodes())
        {
          vertex_id_t nexts = x.get_next_sibling_id();
          if (not pred[nexts])
            {
              nexts = data->m_nodes[nexts].m_next;
              data->m_nodes[x.id()].m_next_sibling = nexts;
            }
          if (not pred[x.id()]) {
            data->m_nodes[x.get_prev_node_id()].m_next = x.get_next_node_id();
            data->m_nodes[x.get_next_node_id()].m_prev = x.get_prev_node_id();
          }

          // Take care of the sentinel
          assert(x.get_parent_id() == tree.npos() or pred[x.get_parent_id()]);
          assert(x.get_next_node_id() == tree.npos() or pred[x.get_next_node_id()]);
          assert(x.get_next_sibling_id() == tree.npos() or pred[x.get_next_sibling_id()]);
        }

      // reaffect point to nodes
      {
        mln_foreach (auto& v, data->m_pmap.values())
          if (not pred[v])
            v = data->m_nodes[v].m_parent;
      }

      mln_exiting();
    }




    template <class P, class Amap, class PredicateMap, class ValueMap, class OutputImage>
    void
    filter_direct_and_reconstruct(const component_tree<P, Amap>& tree,
                                  const PredicateMap& predicate_map,
                                  const ValueMap& value_map,
                                  OutputImage&& out)
    {
      mln_entering("mln::morpho::filter_direct_and_reconstruct");

      typedef component_tree<P, Amap> tree_t;
      typedef typename tree_t::vertex_id_t vertex_id_t;

      //ValueMap vmap = value_map;
      property_map<component_tree<P, Amap>, typename ValueMap::value_type> vmap(tree);

      mln_foreach(auto x, tree.nodes())
        if (not predicate_map[x])
          vmap[x] = vmap[x.parent()];
        else
          vmap[x] = value_map[x];

      mln_foreach(auto px, out.pixels())
        {
          auto x = tree.get_node_at(px.index());
          px.val() = vmap[x];
        }

      mln_exiting();
    }

    template <class P, class Amap, class PredicateMap, class ValueMap>
    property_map<component_tree<P, Amap>, typename ValueMap::value_type>
    filter_direct(const component_tree<P, Amap>& tree,
                  const PredicateMap& predicate_map,
                  const ValueMap& value_map)
    {
      mln_entering("mln::morpho::filter_direct");

      typedef component_tree<P, Amap> tree_t;
      typedef typename tree_t::vertex_id_t vertex_id_t;

      property_map<component_tree<P, Amap>, typename ValueMap::value_type> vmap(tree);

      mln_foreach(auto x, tree.nodes())
        if (not predicate_map[x])
          vmap[x] = vmap[x.parent()];
        else
          vmap[x] = value_map[x];

      mln_exiting();
      return vmap;
    }


  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_FILTERING_HPP
