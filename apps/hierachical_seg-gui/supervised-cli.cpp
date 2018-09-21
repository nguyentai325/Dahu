#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/colors/literal.hpp>

#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/accu/accumulators/mean.hpp>

#include <apps/tos/Kinterpolate.hpp>
#include "myheap.hpp"


using namespace mln;

typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;

property_map<tree_t, float> run_djikstra(const tree_t& tree,
                                         property_map<tree_t, uint8>& colormap,
                                         const property_map<tree_t, rgb<int> >& vmap,
                                         int color)
{
  property_map<tree_t, float> distancemap(tree, value_traits<float>::max());
  property_map<tree_t, int>   posmap(tree, -1);

  myheap heap(tree, distancemap, posmap);
  mln_foreach(auto node, tree.nodes())
    if (colormap[node] == (uint8)-1)
      colormap[node] = 0;
    else if (colormap[node] == color) {
      distancemap[node] = 0;
      heap.push(node);
    }


  {
    mln_entering("Running djikstra");
    tree_t::node_type p;
    morpho::tree_neighb_t nbh;
    mln_iter(q, nbh(p));
    while (not heap.empty())
      {
        p = heap.pop();
        assert(posmap[p] == -1);
        mln_forall(q)
        {
          float d = l2norm(vmap[p] - vmap[*q]) + distancemap[p];
          if (d < distancemap[*q])
            {
              distancemap[*q] = d;
              if (posmap[*q] == -1) { // not in queue, insert
                heap.push(*q);
              } else {
                heap.update(*q);
              }
            }
        }
      }
    mln_exiting();
  }

  return distancemap;
}


image2d<rgb8>
segmentation_(const tree_t& tree,
              const image2d<rgb8>& ima_,
              const image2d<rgb8>& markers__)
{
  image2d<uint8> markers_ = transform(markers__,
                                      [](const rgb8& v) -> uint8 {
                                        if (v == colors::literal::red) return 1;
                                        else if (v == colors::literal::blue) return 2;
                                        else return 0;
                                      });

  image2d<rgb8>  ima = Kadjust_to(ima_, tree._get_data()->m_pmap.domain());
  image2d<uint8> marker = Kadjust_to(markers_, tree._get_data()->m_pmap.domain(), "zero");

  auto vmap = morpho::vaccumulate_proper(tree, ima, accu::features::mean<>());

  property_map<tree_t, uint8> colormap(tree, 0);


  property_map<tree_t, float> distancemap(tree, value_traits<float>::max());
  property_map<tree_t, int>   posmap(tree, -1);

  mln_foreach(auto px, marker.pixels())
    {
      tree_t::node_type node = tree.get_node_at(px.index());
      if (colormap[node] == 0)
        colormap[node] = px.val();
      else if (colormap[node] != px.val())
        colormap[node] = -1;
    }


  myheap heap(tree, distancemap, posmap);
  mln_foreach(auto node, tree.nodes())
    if (colormap[node] == (uint8)-1)
      colormap[node] = 0;
    else if (colormap[node] != 0) {
      distancemap[node] = 0;
      heap.push(node);
    }

  // Run djisktra
  {
    mln_entering("Running djikstra");
    tree_t::node_type p;
    morpho::tree_neighb_t nbh;
    mln_iter(q, nbh(p));
    while (not heap.empty())
      {
        p = heap.pop();
        assert(posmap[p] == -1);
        mln_forall(q)
        {
          float d = l2norm(vmap[p] - vmap[*q]) + distancemap[p];
          if (d < distancemap[*q])
            {
              colormap[*q] = colormap[p];
              distancemap[*q] = d;
              if (posmap[*q] == -1) { // not in queue, insert
                heap.push(*q);
              } else {
                heap.update(*q);
              }
          }
        }
      }
    mln_exiting();
  }

  image2d<float> dist;
  image2d<uint8> color;
  resize(dist, ima);
  resize(color, ima);

  morpho::reconstruction(tree, colormap, color);
  morpho::reconstruction(tree, distancemap, dist);
  io::imsave(color, "/tmp/colormap.tiff");
  io::imsave(dist, "/tmp/distancemap.tiff");

  auto rec = transform(imzip(color, ima), [](const std::tuple<uint8, rgb8>& v) {
      return std::get<0>(v) == 2 ? std::get<1>(v) : rgb8(0);
    });

  return Kadjust_to(rec, ima_.domain());
}


image2d<bool>
segmentation(const tree_t& tree,
             const image2d<rgb8>& ima_,
             const image2d<uint8>& markers__)
{
  image2d<uint8> markers_ = transform(markers__,
                                      [](const uint8& v) -> uint8 {
                                        if (v == 0) return 0; // No tag
                                        else if (v == 255) return 2; // Object
                                        else return 1; // Bg
                                      });

  image2d<rgb8>  ima = Kadjust_to(ima_, tree._get_data()->m_pmap.domain());
  image2d<uint8> marker = Kadjust_to(markers_, tree._get_data()->m_pmap.domain(), "zero");

  auto vmap = morpho::vaccumulate_proper(tree, ima, accu::features::mean<>());

  property_map<tree_t, uint8> colormap(tree, 0);
  mln_foreach(auto px, marker.pixels())
    {
      tree_t::node_type node = tree.get_node_at(px.index());
      if (colormap[node] == 0)
        colormap[node] = px.val();
      else if (colormap[node] != px.val())
        colormap[node] = -1;
    }

  property_map<tree_t, float> dist1 = run_djikstra(tree, colormap, vmap, 1);
  property_map<tree_t, float> dist2 = run_djikstra(tree, colormap, vmap, 2);

  auto color_map = make_functional_property_map<tree_t::node_type>
    ([&dist1,&dist2](const tree_t::node_type& x) -> bool{
      return dist2[x] <= dist1[x];
    });

  auto proba_map = make_functional_property_map<tree_t::node_type>
    ([&dist1,&dist2](const tree_t::node_type& x) {
      return dist1[x] / (dist1[x] + dist2[x]);
    });

  image2d<bool> color;
  resize(color, ima);
  morpho::reconstruction(tree, color_map, color);

  image2d<float> dist;
  resize(dist, ima);
  morpho::reconstruction(tree, proba_map, dist);
  io::imsave(dist, "/tmp/distancemap.tiff");

  return Kadjust_to(color, ima_.domain());
}



int main(int argc, char** argv)
{
  if (argc < 5) {
    std::cout << "Usage: " << argv[0] << "tree input.ppm markers.pgm output.pbm" << std::endl;
    std::exit(1);
  }

  tree_t tree;
  {
    std::ifstream fs(argv[1]);
    morpho::load(fs, tree);
  }

  image2d<rgb8> ima;
  io::imread(argv[2], ima);

  image2d<uint8> markers;
  io::imread(argv[3], markers);

  image2d<bool> res = segmentation(tree, ima, markers);

  io::imsave(res, argv[4]);
}
