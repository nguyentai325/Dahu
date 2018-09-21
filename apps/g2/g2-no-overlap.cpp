#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/always.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/data/stretch.hpp>

#include <mln/morpho/tos/ctos.hpp>
//#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>

#include <mln/accu/accumulators/count.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include <apps/llview/cllview.hpp>
#include <boost/format.hpp>
#include "compute_g2.hpp"
#include "types.hpp"

void usage(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " input.ppm λ α output.tree output.tiff\n"
    "1. Compute the margianal ToS\n"
    "2. Compute the grain filter of size λ\n"
    "3. Compute the height selector (keeps only the llines every α levels )\n"
    "4. Compute G₂\n"
    "5. Keep only non-overlaping shapes (in the R.G.B) order\n"
    ;
  std::exit(1);
}


namespace mln
{
  void mystretch(const image2d< vec<unsigned, NTREE> >& f,
                 const char* name)
  {
    image2d<value_t> out;
    resize(out, f);

    for (int k = 0; k < NTREE; ++k)
      data::stretch_to(channel(f, k), channel(out, k));

    io::imsave(out, name);
  }

}


int main(int argc, char** argv)
{
  if (argc < 6)
    usage(argv);


  using namespace mln;


  image2d<value_t> ima;
  io::imread(argv[1], ima);

  image2d<value_t> f = addborder(ima, lexicographicalorder_less<value_t>());
  image2d<value_t> F = immerse_k1(f);

  typedef value_t::value_type V;

  int ordering[] = {3,0};
  unsigned lambda = std::atoi(argv[2]);
  unsigned alpha = std::atof(argv[3]);
  unsigned alpha2 = std::atof(argv[4]);
  unsigned alphas[] = {alpha, alpha, alpha, alpha2};


  /// Compute the marginal ToS
  tree_t t[NTREE];
  property_map<tree_t, V> val[NTREE];

  for (int i = 0; i < NTREE; ++i)
    {
      auto g = imtransform(f, [i](value_t x) { return x[i]; });
      t[i] = morpho::cToS(g, c4);

      // Grain filter
      {
        auto area = morpho::accumulate(t[i], accu::features::count<> ());
        auto pred = make_functional_property_map<tree_t::vertex_id_t>
          ( [&area, lambda] (tree_t::vertex_id_t x) { return area[x] > lambda; });

        morpho::filter_direct_inplace(t[i], pred);
      }

      // Level lines simplification
      {
        property_map<tree_t, bool> keep(t[i], false);
        val[i] = make_attribute_map_from_image(t[i], channel(F,i));
        auto& v = val[i];

        keep[t[i].get_root()] = true;
        mln_foreach(auto x, t[i].nodes_without_root())
          {
            if ( mln::abs(v[x.parent()] - v[x]) > alphas[i])
              keep[x] = true;
            else
              v[x] = v[x.parent()];
          }

        filter_direct_inplace(t[i], keep);
      }
    }

  // Filter and compress trees
  {
    image2d<V> rec = imchvalue<V>(F);
    for (int i = 0; i < NTREE; ++i)
      {
        morpho::reconstruction(t[i], val[i], rec);
        io::imsave(rec, (boost::format("comp%i.tiff") % i).str());

        t[i]._reorder_pset();
        t[i].shrink_to_fit();
        val[i] = make_attribute_map_from_image(t[i], channel(F,i));
      }
  }

  /// Compute the graph
  typedef Graph<NTREE> MyGraph;
  typedef graph_content<NTREE> my_graph_content;

  MyGraph g2;
  std::array<property_map<tree_t, MyGraph::vertex_descriptor>, NTREE> tlink;
  std::tie(g2, tlink) = compute_g2<NTREE>(t);


  /// Algo to select non-overlapping RGB

  auto glink = boost::get(&my_graph_content::tlinks, g2);
  auto d = boost::get(&my_graph_content::depth, g2);


  typedef mln::vec<unsigned, NTREE> depth_vec_t;

  image2d< mln::vec<unsigned, NTREE> > imdepth;
  resize(imdepth, F);

  image2d<rgb8> selection;
  image2d<rgb8> rec;
  resize(selection, F);
  resize(rec, F);

  // Debug
  {
    for (int k = 0; k < NTREE; ++k)
      {
        auto dmap = make_functional_property_map<tree_t::node_type>([&d,&tlink,k](tree_t::node_type x) {
            return d[tlink[k][x]][k];
          });
        morpho::reconstruction(t[k], dmap, channel(imdepth,k));
      }
    io::imsave(imcast<value_t>(imdepth), "depth.tiff");
  }

  // 1st step
  {
    int k = ordering[0];
    auto dmap = make_functional_property_map<tree_t::node_type>([&d,&tlink,k](tree_t::node_type x) {
        return d[tlink[k][x]];
      });
    morpho::reconstruction(t[k], dmap, imdepth);

    // Show those llines on selection map
    auto smap = make_functional_property_map<tree_t::node_type>(yes_t ());
    auto vmap = make_functional_property_map<tree_t::node_type>(always_t<uint8> (255));

    copy(llview_val(t[k], dmap, smap, vmap), channel(selection, 0));
    if (k != 3) {
      morpho::filter_direct_and_reconstruct(t[k], smap, val[k], channel(rec, k));
    }

    io::imsave(imcast<value_t>(imdepth), (boost::format("select%i.tiff") % k).str());
  }

  // 2st step
  // Accumulate the min/max elements
  for (int i = 1; i < sizeof(ordering) / sizeof(int); ++i)
    {
      int k = ordering[i];

      property_map<tree_t, bool> pred(t[k], true);
      mln_foreach(auto px, imdepth.pixels())
        {
          tree_t::node_type x = t[k].get_node_at(px.index());
          depth_vec_t v = px.val();
          while (x.id() != tree_t::npos())
            {
              depth_vec_t val = d[tlink[k][x]];
              pred[x] = pred[x] and (vecprod_islessequal(v, val) or
                                     vecprod_islessequal(val, v));
              x = x.parent();
            }
        }


      auto dmap = make_functional_property_map<tree_t::node_type>([&d,&tlink,k](tree_t::node_type x) {
          return d[tlink[k][x]];
        });


      {
        auto vmap = make_functional_property_map<tree_t::node_type>(always_t<uint8> (255));
        if (i < 3)
          copy(llview_val(t[k], dmap, pred, vmap), channel(selection, i));

        if (k != 3) {
          morpho::filter_direct_and_reconstruct(t[k], pred, val[k], channel(rec, k));
        }
      }

      {
        unsigned n = 0;
        mln_foreach (auto x, t[k].nodes()) {
          //std::cout << dmap[x];
          if (pred[x])
            n++;
          // else
          //   std::cout << " K";
          // std::cout << std::endl;
        }
        std::cout << "Selected " << n << " over " << t[k].size() << " from tree " << k << std::endl;
      }

      image2d< mln::vec<unsigned, NTREE> > tmp;
      resize(tmp, F);
      morpho::filter_direct_and_reconstruct(t[k], pred, dmap, tmp);
      imdepth = eval(isup(imdepth, tmp));
      io::imsave(imcast<value_t>(tmp), (boost::format("select%i.tiff") % k).str());
    }

  io::imsave(selection, "selection.tiff");
  io::imsave(rec, "reconstruction.tiff");
}
