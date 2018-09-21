#include <mln/core/image/image2d.hpp>

#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>
#include <mln/morpho/extinction.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <apps/tos/Kinterpolate.hpp>
#include <apps/saliency/attribute_on_contour.hpp>

#include <fstream>
#include "gradient_magnitude.hpp"

int main(int argc, char** argv)
{
  if (argc != 6)
    {
      std::cerr << "Usage:" << argv[0] << " tree input.tiff λ t₁ out.tiff\n"
        "Export a saliency map from this attribute\n"
        " λ:\t Grain filter before anything else (number of 2F)\n"
        " t₁:\t Threshold above which node having an energy greater that t₁ are removed.\n";
      std::exit(1);
    }

  using namespace mln;


  const char* tree_path = argv[1];
  const char* img_path = argv[2];
  unsigned grain = std::atoi(argv[3]);
  float threshold1 = std::atof(argv[4]);
  const char* output_path = argv[5];


  typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;
  tree_t tree;
  {
    std::ifstream f(tree_path, std::ios::binary);
    morpho::load(f, tree);
  }

  typedef uint8 V;

  image2d<V> ima_, ima;
  io::imread(img_path, ima_);

  ima = Kadjust_to(ima_, tree._get_data()->m_pmap.domain());

  if (ima.domain() != tree._get_data()->m_pmap.domain())
    {
      std::cerr << "Domain between image differ.\n"
                << ima.domain() << " vs "
                << tree._get_data()->m_pmap.domain() << std::endl;
      std::exit(1);
    }

  // 1. grain filter
  {
    accu::accumulators::accu_if< accu::accumulators::count<>,
                                 K1::is_face_2_t,
                                 point2d > counter;
    auto areamap = morpho::paccumulate(tree, ima, counter);

    auto pred = make_functional_property_map<tree_t::vertex_id_t>([&areamap, grain](tree_t::vertex_id_t x) {
        return areamap[x] > grain;
      });
    morpho::filter_direct_inplace(tree, pred);
  }

  // 2. Compute energy
  auto energy = compute_gradient_magnitude(tree, ima);

  // 3. Convert to image and filter
  auto ienergy = morpho::make_image(tree, energy);
  {
    mln_foreach(float& v, ienergy.values())
      if (v > threshold1)
        v = threshold1;
  }

  // 4. extinction
  auto extincted = morpho::extinction(ienergy, morpho::tree_neighb_t());

  // 5. reconstruction
  {
    auto& attr = extincted.get_vmap();
    image2d<float> sal = imchvalue<float>(ima).init(0);
    attribute_on_contour(tree, attr, sal);
    io::imsave(sal, output_path);
  }

}
