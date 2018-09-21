#include <mln/io/imread.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/accu/accumulators/accu_if.hpp>
#include <mln/accu/accumulators/count.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include "attributes.hpp"
#include "cMSER.hpp"



int main(int argc, char** argv)
{
  using namespace mln;

  po::options_description desc("MSER Options test");
  desc.add_options()
    ("eps,h", po::value<float>()->required(), "The height for neighbor lookup (e.g. 10).");

  po::variables_map vm = process_cmdline(argc, argv, desc);
  //tree_t tree = preprocess(vm);



  image2d<float> f;
  io::imread(vm["input_path"].as<std::string>(), f, true);

  io::imsave(f, "test.ppm");

  // image2d<float> F = Kadjust_to(f, tree._get_data()->m_pmap.domain());

  // auto vmap = morpho::make_attribute_map_from_image(tree, F);

  // accu::accumulators::accu_if< accu::accumulators::count<>,
  //                              K1::is_face_2_t,
  //                              point2d > counter;
  // auto amap = morpho::paccumulate(tree, F, counter);

  // auto amser = compute_MSER(tree, vmap, amap, vm["eps"].as<float>(), MSER_NORM);

  // auto smap = postprocess(vm, tree, amser);

  // std::cout << "test" << std::endl;

  // const char* names[] = {"MSER", "Extinction"};
  // std::function<float(tree_t::node_type)> attrs[] = { _as_fun(amser), _as_fun(smap) };

  // export_(vm, tree, smap, names, attrs, 2);
}
