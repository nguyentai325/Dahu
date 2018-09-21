#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/accu/accumulators/count.hpp>
#include <mln/morpho/maxtree/maxtree.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>

#include <mln/io/imprint.hpp>
#include <iostream>

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 5)
    {
      std::cerr << "Usage: " << argv[0] << " input(gray) connectivity lambda output(gray)" << std::endl
		<< "connectivity: 4 / 8" << std::endl;
      std::exit(1);
    }

  image2d<uint8> ima;
  io::imread(argv[1], ima);

  typedef image2d<uint8>::size_type size_type;
  morpho::component_tree< size_type, image2d<size_type> > tree;

  if (argv[2][0] == '4')
    tree = morpho::maxtree_indexes(ima, c4, std::greater<uint8>() );
  else
    tree = morpho::maxtree_indexes(ima, c8, std::greater<uint8>());

  auto areamap = morpho::paccumulate(tree, ima, accu::features::count<> ());

  unsigned lambda = std::atoi(argv[3]);
  auto criterion = make_functional_property_map<size_type>
    ([&areamap, lambda](size_type x) {
      return areamap[x] >= lambda;
    });

  morpho::filter_direct_inplace(tree, criterion);

  auto valuemap = morpho::make_attribute_map_from_image(tree, ima);

  image2d<uint8> out;
  resize(out, ima);
  morpho::reconstruction(tree, valuemap, out);

  io::imsave(out, argv[4]);
}
