#include <iostream>
#include <mln/core/image/image2d.hpp>
#include <mln/core/colors.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/morpho/component_tree/component_tree.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/accu/accumulators/mean.hpp>
#include <mln/accu/accumulators/accu_if.hpp>
#include <mln/accu/accumulators/accu_transform.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>

int main(int argc, char** argv)
{
  if (argc < 4)
    {
      std::cerr << "Usage: " << argv[0] << " input.tree img(color|gray)  out.tiff \n"
                << "Reconstruct from a tree. The image may in K0|K1 with(out) border.\n"
        ;
      std::exit(1);
    }

  using namespace mln;


  morpho::component_tree<unsigned, image2d<unsigned> > tree;
  image2d<rgb8> ima;

  {
    std::ifstream fs(argv[1], std::ios::binary);
    morpho::load(fs, tree);
  }

  io::imread(argv[2], ima);

  image2d<rgb8> f = Kadjust_to(ima, tree._get_data()->m_pmap.domain());

  typedef accu::accumulators::mean<rgb8> ACCU;
  typedef image2d<rgb8> I;

  auto is_face_2 = [](const mln_cpixel(I)& px) { return K1::is_face_2(px.point()); };
  auto val_of_pix = [](const mln_cpixel(I)& px) { return px.val(); };

  typedef accu::accumulators::accu_transform<ACCU, decltype(val_of_pix), mln_cpixel(I)> ACCU_1;
  typedef accu::accumulators::accu_if<ACCU_1, decltype(is_face_2), mln_cpixel(I)> ACCU_2;

  //auto vmap = morpho::pixaccumulate_proper(tree, f, ACCU_2(ACCU_1(ACCU(), val_of_pix),is_face_2));
  auto vmap = morpho::vaccumulate_proper(tree, f, accu::features::mean<>());

  image2d<rgb8> out;
  resize(out, f);
  morpho::reconstruction(tree, vmap, out);

  image2d<rgb8> final = Kadjust_to(out, ima.domain());
  io::imsave(final, argv[3]);
}




