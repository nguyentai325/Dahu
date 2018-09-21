#include "attributes.hpp"
#include <apps/tos/croutines.hpp>
#include <mln/morpho/extinction.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/filtering.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/io/imsave.hpp>

namespace mln
{


  po::variables_map
  process_cmdline(int argc, char** argv, po::options_description desc_, const char* usage)
  {
    po::options_description hidden("Allowed options");
    hidden.add_options()
      ("tree_path",  po::value<std::string>()->required(), "Input tree")
      ("input_path", po::value<std::string>()->required(), "Input file (gray level)")
    ;

    po::options_description desc("General options");
    desc.add_options()
      ("help", "Help message")
      ("grain,g", po::value<int>(), "Perform a grain filter before anything else");

    po::options_description desc2("Extinction Options");
    desc2.add_options()
      ("clip", po::value<float>(), "Clip the energy map to [0, t₁]");

    po::options_description desc3("Output Options");
    desc3.add_options()
      ("output-tree", po::value<std::string>(), "Path to save the tree.")
      ("output-saliency", po::value<std::string>(), "Path to save the saliency map")
      ("output-csv", po::value<std::string>(), "Path to the csv file to save the attributes")
      ("output-simp-tree", po::value<std::string>(), "Path to save the simplified tree.")
      ("output-simp-depth", po::value<std::string>(), "Path to save the simplified tree depth image.")
      ;

    po::positional_options_description pd;
    pd.add("tree_path", 1)
      .add("input_path", 1)
      ;

    po::options_description all("Allowed options");
    po::options_description visible("Allowed options");
    all.add(desc_).add(hidden).add(desc).add(desc2).add(desc3);
    visible.add(desc_).add(desc).add(desc2).add(desc3);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                .options(all)
                .positional(pd).run(), vm);

      po::notify(vm);
    } catch (...) {
      std::cerr << "Usage: " << argv[0] << " input.tree input.tiff [options]\n"
                << usage
                << visible;
      std::exit(1);
    }

    if (vm.count("help"))
      {
        std::cerr << "Usage: " << argv[0] << " input.tree input.tiff [options]\n"
                  << usage
                  << visible;
        std::exit(1);
    }

    return vm;
  }


  tree_t
  preprocess(const po::variables_map& vm)
  {
    tree_t tree;
    morpho::load(vm["tree_path"].as<std::string>(), tree);

    if (vm.count("grain"))
      grain_filter_inplace(tree, vm["grain"].as<int>());

    if (vm.count("output-tree"))
      morpho::save(tree, vm["output-tree"].as<std::string>());

    return tree;
  }


  property_map<tree_t, float>
  postprocess(const po::variables_map& vm,
              const tree_t& tree,
              const property_map<tree_t, float>& energy)
  {
    // convert to image and filter
    auto ienergy = morpho::make_image(tree, energy);
    if (vm.count("clip"))
    {
      float threshold1 = vm["clip"].as<float>();
      mln_foreach(float& v, ienergy.values())
        if (v > threshold1)
          v = threshold1;
    }

    auto extincted = morpho::extinction(ienergy, morpho::tree_neighb_t()// , std::greater<float>()
      );

    property_map<tree_t, float> vmap = std::move(extincted.get_vmap());
    return vmap;
  }

  void export_(const po::variables_map& vm,
               tree_t& tree,
               const property_map<tree_t, float>& saliency_map,
               const char* names[],
               std::function<float(tree_t::node_type)> vmaps[],
               int sz)
  {
    mln_entering("export_");

    {
      int nnodes = 0;
      mln_foreach(auto x, tree.nodes())
        if (saliency_map[x] > 0)
          nnodes++;

      std::cout << "Remaining nodes: " << nnodes << std::endl;
    }


    if (vm.count("output-saliency"))
      {
        image2d<float> sal = set_value_on_contour(tree, saliency_map);
        io::imsave(sal, vm["output-saliency"].as<std::string>());
      }

    if (vm.count("output-csv"))
      {
        std::ofstream f(vm["output-csv"].as<std::string>());
        for (int i = 0; i < sz; ++i)
          f << names[i] << ((i+1 < sz) ? ',' : '\n');

        mln_foreach(auto x, tree.nodes())
          for (int i = 0; i < sz; ++i)
            f << vmaps[i](x) << ((i+1 < sz) ? ',' : '\n');
      }

    if (vm.count("output-simp-tree") or
        vm.count("output-simp-depth") or
        vm.count("export-rec"))
      {
        auto binmap = make_functional_property_map<tree_t::vertex_id_t>
          ([&saliency_map](tree_t::vertex_id_t x) -> bool {
            return saliency_map[x] > 0;
          });


        morpho::filter_direct_inplace(tree, binmap);

        if (vm.count("output-simp-tree"))
          morpho::save(tree, vm["output-simp-tree"].as<std::string>());

        if (vm.count("output-simp-depth"))
          {
            auto depth = morpho::compute_depth(tree);
            image2d<uint16> d;
            d.resize(tree._get_data()->m_pmap.domain());
            morpho::reconstruction(tree, depth, d);
            io::imsave(d, vm["output-simp-depth"].as<std::string>());
          }
      }
    mln_exiting();
  }



}
