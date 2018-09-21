#ifndef MLN_MORPHO_ALPHATREE_ALPHATREE_HPP
# define MLN_MORPHO_ALPHATREE_ALPHATREE_HPP

# include <mln/morpho/component_tree/component_tree.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>

namespace mln
{

  namespace morpho
  {

    template <class I, class N, class Distance = functional::l2dist_t<mln_value(I)> >
    std::pair<
      morpho::component_tree<typename I::size_type, mln_ch_value(I, typename I::size_type)>,
      property_map<
        morpho::component_tree<typename I::size_type, mln_ch_value(I, typename I::size_type)>,
        typename std::result_of<Distance(mln_value(I), mln_value(I))>::type
        >
      >
    alphatree_indexes(const Image<I>& f, const Neighborhood<N>& nbh, Distance d = Distance ());


    /*******************************/
    /***   Implementation        ***/
    /*******************************/

    namespace internal
    {
      inline
      unsigned
      _zfindroot(std::vector<unsigned>& par, unsigned i)
      {
        if (i == par[i])
          return i;
        else
          return par[i] = _zfindroot(par, par[i]);
      }
    }


    template <class I, class N, class Distance>
    std::pair<
      morpho::component_tree<typename I::size_type, mln_ch_value(I, typename I::size_type)>,
      property_map<
        morpho::component_tree<typename I::size_type, mln_ch_value(I, typename I::size_type)>,
        typename std::result_of<Distance(mln_value(I), mln_value(I))>::type
        >
      >
    alphatree_indexes(const Image<I>& f_, const Neighborhood<N>& nbh_, Distance distance)
    {
      mln_entering("mln::morpho::alphatree_indexes");

      const I& f = exact(f_);
      const N& nbh = exact(nbh_);

      typedef typename I::size_type index_t;
      typedef morpho::component_tree<index_t, mln_ch_value(I, index_t)> tree_t;
      typedef typename std::result_of<Distance(mln_value(I), mln_value(I))>::type distance_t;

      std::vector< std::tuple<index_t, index_t, distance_t> > edges;

      // 1st step: store and sort edges
      {
        mln_pixter(px, f);
        mln_iter(nx, nbh(px));

        edges.reserve(f.domain().size());

        mln_forall(px)
        {
          mln_forall(nx)
            if (f.domain().has(nx->point()) and px->index() < nx->index())
              edges.emplace_back(px->index(), nx->index(), distance(px->val(), nx->val()));
        }



        std::sort(edges.begin(), edges.end(),
                  [distance,&f](std::tuple<index_t, index_t, distance_t> v1,
                                std::tuple<index_t, index_t, distance_t> v2) {
                    return std::get<2>(v1) < std::get<2>(v2);
                  });
      }

      // 2nd step: make the tree
      tree_t                            tree;
      property_map<tree_t, distance_t>  vmap;
      {
        mln_ch_value(I, unsigned) zpar;
        std::vector<unsigned>     par;
        std::vector<unsigned>     vzpar;
        std::vector<distance_t>   dist;
        std::size_t sz = f.domain().size();

        par.resize(sz);
        vzpar.resize(sz);
        dist.resize(sz);
        resize(zpar, f);

        // 2.1 Make set for each pixel
        unsigned k = 0;
        mln_foreach(auto px, zpar.pixels())
          {
            px.val() = k;
            par[k] = k;
            vzpar[k] = k;
            dist[k] = 0;
            ++k;
          }

        // 2.2 union find
        for (auto e : edges)
          {
            index_t x,y;
            distance_t d;
            std::tie(x,y,d) = e;
            unsigned rx = (zpar[x] = internal::_zfindroot(vzpar, zpar[x]));
            unsigned ry = (zpar[y] = internal::_zfindroot(vzpar, zpar[y]));

            // std::cout << "Process: " << x << " " << y << " dist: " << d << "\n"
            //           << "   " << (int)f[x] << "  " << (int)f[y] << "\n";

            if (rx != ry) // Merge set
              {
                // std::cout << "   Merging " << rx << " and " << ry << "\n";
                par.push_back(k);
                vzpar.push_back(k);
                dist.push_back(d);
                zpar[x] = vzpar[rx] = par[rx] = k;
                zpar[y] = vzpar[ry] = par[ry] = k;
                ++k;
              }
            // else
            //   {
            //     std::cout << "   Nothing to merge\n";
            //   }
          }

        // 2.3 Canonicalization
        unsigned nnodes = 0;
        {
          for (int p = k-1; p >= 0; --p)
            {
              unsigned q = par[p];
              if (dist[q] != dist[p] or q == par[q]) // root case
                ++nnodes;

              if (dist[par[q]] == dist[q])
                par[p] = par[q];
            }
          unsigned p = 0;
          mln_foreach(auto px, zpar.pixels())
            {
              unsigned q = par[p];
              if (dist[q] == dist[p]) {
                px.val() = q;
              } else {
                px.val() = p;
              }
              ++p;
            }
        }

        // 2.4 Build the tree from the parent image
        auto data = tree._get_data();
        auto& nodes = data->m_nodes;
        data->m_nodes.resize(nnodes + 1);
        data->m_S.resize(sz);
        resize(data->m_pmap, f);

        vmap = property_map<tree_t, distance_t> (tree);

        k = 1;
        // Handle sentinel and root
        data->m_nodes[1] = internal::component_tree_node {
          tree_t::npos(), // parent
          tree_t::npos(), // prev
          tree_t::npos(), // next
          tree_t::npos(), // nexts
          0,            // point index (Start of buffer)
        };
        data->m_nodes[tree_t::npos()] = internal::component_tree_node {
          tree_t::npos(), // parent
          1,              // prev = last_node = root
          tree_t::npos(), // next = itself
          tree_t::npos(), // next sibling -> itself
          (unsigned)sz    // point index (End of Buffer)
        };

        vmap[1] = dist.back();  // set the root value
        par.back() = 1;         // par points to the new index
        ++k;

        for (int i = par.size() - 2; i >= 0; --i)
          {
            if (dist[i] == dist[par[i]])
              {
                par[i] = par[par[i]];
                continue;
              }

            // Insert the node in child of par[i]
            unsigned q = par[par[i]];
            unsigned next = nodes[q].m_next;
            nodes[k].m_parent = q; // set parent
            // insert in the dble linked list
            nodes[next].m_prev = k;
            nodes[q].m_next = k;
            nodes[k].m_prev = q;
            nodes[k].m_next = next;
            // insert in the skip list
            nodes[k].m_next_sibling = next;

            // update par
            par[i] = k;
            vmap[k] = dist[i];
            ++k;
          }

        // Insert in S and set pmap relation
        {
          unsigned p = 0;
          unsigned i = sz-1;
          mln_foreach(auto px, data->m_pmap.pixels())
            {
              data->m_S[i] = px.index();
              px.val() = par[p];
              nodes[par[p]].m_point_index = i;
              ++p;
              --i;
            }
        }
      }

      mln_exiting();
      return  { std::move(tree.get_subtree(1)), std::move(vmap) };
    }

  }

}

#endif // ! MLN_MORPHO_ALPHATREE_ALPHATREE_HPP
