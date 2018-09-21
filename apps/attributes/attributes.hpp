#ifndef APPS_ATTRIBUTE_ATTRIBUTES_HPP
# define APPS_ATTRIBUTE_ATTRIBUTES_HPP

# include <boost/program_options.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>

namespace po = boost::program_options;
typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

namespace mln
{

  po::variables_map
  process_cmdline(int argc, char** argv, po::options_description desc, const char* usage = "");

  tree_t
  preprocess(const po::variables_map& vm);

  property_map<tree_t, float>
  postprocess(const po::variables_map& vm,
              const tree_t& tree,
              const property_map<tree_t, float>& energy);



  void export_(const po::variables_map& vm,
               tree_t& tree,
               const property_map<tree_t, float>& saliency_map,
               const char* names[],
               std::function<float(tree_t::node_type)> vmaps[],
               int sz);


  template <class Amap>
  std::function<float(tree_t::node_type)> _as_fun(const Amap& w) {
    return [&w](tree_t::node_type x) -> float { return w[x]; };
  }


}

#endif // ! APPS_ATTRIBUTE_ATTRIBUTES_HPP
