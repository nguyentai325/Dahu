#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/win2d.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/literal.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/morpho/alphatree/alphatree.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>
#include <mln/morpho/component_tree/cuts.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/accu/accumulators/infsup.hpp>
#include <mln/accu/accumulators/mean.hpp>
#include <mln/transform/chamfer_distance_transform.hpp>
#include <mln/morpho/algebraic_filter.hpp>

int main(int argc, char** argv)
{
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] << " input.ppm BV grain #markers markers.pgm" << std::endl;
    std::exit(1);
  }

  using namespace mln;

  typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;
  typedef rgb8 V;

  std::string outname = argv[5];
  std::string stem = outname.substr(0, outname.rfind('.'));

  image2d<V> f;
  io::imread(argv[1], f);



  tree_t tree;
  property_map<tree_t, double> vmap;

  std::tie(tree, vmap) = morpho::alphatree_indexes(f, c4);

  auto infsup = morpho::vaccumulate(tree, f, accu::accumulators::infsup<V> ());
  auto bvmap = make_functional_property_map<tree_t::node_type>
    ([&infsup](const tree_t::node_type& x) { return l2dist(infsup[x].first, infsup[x].second); });
  auto bvmap2 = make_functional_property_map<unsigned>
    ([&infsup](unsigned x) { return l2dist(infsup[x].first, infsup[x].second); });

  // Filter the alpha tree
  int maxbv = std::atoi(argv[2]);
  morpho::cut_inplace(tree,
                      make_functional_property_map<tree_t::node_type>
                      ([bvmap, maxbv](tree_t::node_type x) {
                        return bvmap[x] < maxbv;
                      }));


  

  image2d<V> out;
  resize(out, f);

  auto avmap = morpho::vaccumulate(tree, f, accu::features::mean<>());
  auto mc = morpho::paccumulate(tree, f, accu::accumulators::mean<point2df>());
  auto areamap = morpho::accumulate(tree, accu::features::count<>());

  morpho::reconstruction(tree, avmap, out);

  std::cout << "See: /tmp/alpha.tiff" << std::endl;
  io::imsave(out, "/tmp/alpha.tiff");

  // Get the seeds.
  unsigned area_threshold = std::atoi(argv[3]);
  property_map<tree_t, uint8> active(tree, true);
  std::vector<unsigned> seeds;

  // Backward
  mln_reverse_foreach(auto x, tree.nodes())
    {
      if (not active[x])
        active[x.parent()] = false;
      else if (areamap[x] > area_threshold)
        {
          active[x.parent()] = false;
          seeds.push_back(x.id());
        }
      else
        {
          active[x] = false;
        }
    }

  auto energy = [&areamap, &bvmap2](unsigned x) -> int { return areamap[x]; }; // - 50 * bvmap2[x]; };
  std::sort(seeds.begin(), seeds.end(), [&energy](unsigned x, unsigned y) { return energy(x) > energy(y); });


  // Only retain the # first seeds.
  int NMarkers = std::min<int>(std::atoi(argv[4]), seeds.size());

  for (int i = 0; i < seeds.size(); ++i)
    active[seeds[i]] = (i < NMarkers) ? i+1 : 0;

  // Forward
  mln_foreach(auto x, tree.nodes())
    if (active[x.parent()])
      active[x] = active[x.parent()]; // label propagation


  // Output the seeds
  image2d<uint8> ske, imlabel;
  resize(ske, f).init(0);
  resize(imlabel, f).init(0);

  {
    // Distance transform
    mln_foreach(auto px, ske.pixels())
      {
        tree_t::node_type x = tree.get_node_at(px.index());
        if (active[x])
          px.val() = 1;
      }

    // Separate the regions
    {
      mln_pixter(px, ske);
      mln_iter(q, c4(px));
      mln_forall(px) {
        unsigned id1 = active[tree.get_node_at(px->index())];
        mln_forall(q) {
          unsigned id2 = active[tree.get_node_at(q->index())];
          if (id1 != id2 and id2 != 0) {
            px->val() = 0;
            break;
          }
        }
      }
    }

    ske = morpho::area_closing(ske, c8, 20);

    transform::chamfer_distance_transform(ske, c4, ske);

    std::cout << "See: " << stem << ".ske.pgm\n";
    io::imsave(ske, stem + ".ske.pgm");


    // Mark the pixels which are far enough from the borders
    auto avg_dist = morpho::vaccumulate(tree, ske, accu::features::mean<float>());
    mln_pixter(px1, px2, ske, imlabel);
    mln_forall(px1, px2)
      {
        tree_t::node_type x = tree.get_node_at(px1->index());
        if (active[x] and px1->val() > avg_dist[x])
          px2->val() = active[x];
      }
  }


  {
    rect2d mywin = make_rectangle2d(11, 11);

    point2d p;
    mln_iter(q, mywin(p));

    for (int i = 0; i < NMarkers; ++i)
      {
        unsigned id = seeds[i];
        std::cout << "Seed #" << i << " : " << energy(id) << std::endl;
        p = mc[id];
        mln_forall(q) {
          if (out.domain().has(*q)) {
            unsigned node_id = tree.get_node_at(f.index_of_point(*q)).id();
            if (node_id == id) {
              out(*q) = colors::literal::blue;
              //imlabel(*q) = i+1;
            }
          }
        }
      }
  }

  // Clean Remove small CC
  imlabel = morpho::area_opening(imlabel, c8, 100);

  io::imsave(imlabel, argv[5]);
  std::cout << "See: " << stem << ".markers.tiff\n";
  io::imsave(out, stem + ".markers.tiff");
}
