#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/transform.hpp>

#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/maxtree/maxtree.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>

using namespace mln;
typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;

///
///
/// \param tree The component tree of the saliency map
/// \param markers The image with markers
image2d<bool>
segmentation(const tree_t& tree,
             const image2d<uint8>& markers__)
{
  enum colortag { BLANC = 0, ROUGE = 1, NOIR = 2, NONE = 3 };

  auto markers = transform(markers__,
                           [](const uint8& v) -> uint8 {
                             if (v == 0) return BLANC; // No tag
                             else if (v == 255) return ROUGE; // Object
                             else return NOIR; // Bg
                           });

  property_map<tree_t, uint8> tags(tree, NONE);
  mln_foreach(auto px, markers.pixels())
    {
      tree_t::node_type x = tree.get_node_at(px.index());
      colortag v = (colortag)px.val();
      if (K1::is_face_2(px.point()) and v != BLANC)
        {
          if (tags[x] == NONE)
            tags[x] = v;
          else if (tags[x] != v)
            tags[x] = BLANC;
        }
    }

  tags[tree.get_root()] = NOIR; // Root is background

  // Propagate up
  mln_reverse_foreach(auto x, tree.nodes()) {
    if (tags[x] != BLANC) {
      auto q = x.parent();
      if (tags[q] == BLANC) // Ok propagate color
        tags[q] = tags[x];
      else if (tags[q] != tags[x]) // Either bg thus both red/black => thus black
        tags[q] = NOIR;
    }
  }

  // Propagate down
  mln_foreach(auto x, tree.nodes_without_root()) {
    auto q = x.parent();
    mln_assertion(tags[q] != BLANC);
    if (tags[x] == BLANC)
      tags[x] = tags[q];
  }

  //
  image2d<bool> out;
  resize(out, markers__);

  auto vmap = make_functional_property_map<tree_t::node_type>
    ( [&tags](const tree_t::node_type& x) { return tags[x] == ROUGE; } );

  morpho::reconstruction(tree, vmap, out);

  return out;
}

int main(int argc, char** argv)
{
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << "saliency(float) markers.pgm output.pbm" << std::endl;
    std::exit(1);
  }

  image2d<float> sal;
  image2d<uint8> markers;
  image2d<bool> out;
  tree_t tree;



  io::imread(argv[1], sal);
  io::imread(argv[2], markers);
  box2d dom = markers.domain();

  markers =  Kadjust_to(markers, sal.domain());
  tree = morpho::mintree_indexes(sal, c4);
  out = segmentation(tree, markers);
  out = Kadjust_to(out, dom);

  io::imsave(out, argv[3]);
}
