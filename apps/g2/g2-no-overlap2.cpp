#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/always.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/data/stretch.hpp>

#include <mln/morpho/tos/ctos.hpp>
//#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>

#include <mln/accu/accumulators/count.hpp>
#include <mln/accu/accumulators/mean.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include <apps/llview/cllview.hpp>
#include <apps/attributes/read.hpp>
#include <boost/format.hpp>
#include "compute_g2.hpp"
#include "types.hpp"

void usage(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " input.tiff channel1.tree channel1.csv channel2.tree channel2.csv\n"
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

  image2d<value_t> f = addborder_marginal(ima);
  image2d<value_t> F = interpolate_k1(f);

  image2d<uint8> F2 = eval(channel(F,1));

  typedef value_t::value_type V;

  enum { NTREE = 2 };
  int chan[NTREE] = {3,1};

  /// Compute the marginal ToS
  tree_t t[2];
  property_map<tree_t, V> val[2];
  property_map<tree_t, float> energy[2];
  property_map<tree_t, int> mval[2];

  morpho::load(argv[2], t[0]);
  energy[0] = read_attribute_map<float>(argv[3], "extinction", t[0]);

  morpho::load(argv[4], t[1]);
  energy[1] = read_attribute_map<float>(argv[5], "extinction", t[1]);

  for (int i = 0; i < 2; ++i)
  {
    tree_t::vertex_id_t root_id = t[i].get_root_id();
    auto pred = make_functional_property_map<tree_t::vertex_id_t>
      ([&energy,root_id,i] (tree_t::vertex_id_t x) { return x == root_id or energy[i][x] > 0; });
    morpho::filter_direct_inplace(t[i], pred);
  }
  io::imsave(f, "tmp.tiff");

  // Filter and compress trees
  {
    image2d<V> rec = imchvalue<V>(F);
    for (int i = 0; i < 2; ++i)
    {
      val[i] = make_attribute_map_from_image(t[i], channel(F, chan[i]));
      morpho::reconstruction(t[i], val[i], rec);
      io::imsave(rec, (boost::format("comp%i.tiff") % i).str());
      t[i]._reorder_pset();
      t[i].shrink_to_fit();
      val[i] = make_attribute_map_from_image(t[i], channel(F, chan[i]));
      mval[i] = vaccumulate_proper(t[i], F2, accu::features::mean<> ());
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
  image2d<uint8> out;
  resize(selection, F);
  resize(rec, F);
  resize(out, F);

  auto cvrter_fun = [] (depth_vec_t x) {
    return rgb8{ x[0], (uint8)((NTREE > 1) ? x[1] : 0), (uint8)((NTREE > 2) ? x[2] : 0) };
  };
  auto cvrter = [cvrter_fun](const image2d<depth_vec_t>& f) {
    return imtransform(f, cvrter_fun);
  };

  // Debug
  {
    for (int k = 0; k < NTREE; ++k)
      {
        auto dmap = make_functional_property_map<tree_t::node_type>([&d,&tlink,k](tree_t::node_type x) {
            return d[tlink[k][x]][k];
          });
        morpho::reconstruction(t[k], dmap, channel(imdepth,k));
      }
    io::imsave(cvrter(imdepth), "depth.tiff");
  }

  // 1st step
  {
    int k = 0;
    auto dmap = make_functional_property_map<tree_t::node_type>([&d,&tlink,k](tree_t::node_type x) {
        return d[tlink[k][x]];
      });
    morpho::reconstruction(t[k], dmap, imdepth);

    // Show those llines on selection map
    auto smap = make_functional_property_map<tree_t::node_type>(yes_t ());
    auto vmap = make_functional_property_map<tree_t::node_type>(always_t<uint8> (255));

    copy(llview_val(t[k], dmap, smap, vmap), channel(selection, 0));

    morpho::filter_direct_and_reconstruct(t[k], smap, mval[k], channel(rec, k));
    copy(channel(rec, k), out);

    io::imsave(cvrter(imdepth), (boost::format("select%i.tiff") % k).str());
  }

  // 2st step
  // Accumulate the min/max elements
  for (int i = 1; i < NTREE; ++i)
    {
      int k = i;

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
          morpho::filter_direct_and_reconstruct(t[k], pred, mval[k], channel(rec, k));
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

      auto dom = where(imzip(imdepth, tmp), [](std::tuple<depth_vec_t, depth_vec_t> v) -> bool {
          return vecprod_isless(std::get<0>(v), std::get<1>(v));
        });
      copy(channel(rec,k) | dom, out | dom);

      imdepth = eval(isup(imdepth, tmp));
      io::imsave(cvrter(tmp), (boost::format("select%i.tiff") % k).str());
    }

  io::imsave(selection, "selection.tiff");
  io::imsave(out, "reconstruction.tiff");
}
