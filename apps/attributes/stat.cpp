#include <mln/core/image/image2d.hpp>
#include <mln/core/colors.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>

#include <mln/accu/accumulators/accu_as_it.hpp>
#include <mln/accu/accumulators/variance.hpp>

int main(int argc, char** argv)
{
  if (argc != 4)
    {
      std::cerr << "Usage:" << argv[0] << " tree input.tiff out.csv\n";
      std::exit(1);
    }

  using namespace mln;

  const char* tree_path = argv[1];
  const char* img_path = argv[2];
  const char* output_path = argv[3];


  typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;
  tree_t tree;
  morpho::load(tree_path, tree);

  image2d<float> ima;
  io::imread(img_path, ima, true);
  ima = Kadjust_to(ima, tree._get_data()->m_pmap.domain());

  auto attr = accu::accumulators::accu_as_it< accu::accumulators::variance<float> >();
  auto res = morpho::vaccumulate(tree, ima, attr);
  auto res2 = morpho::vaccumulate_proper(tree, ima, attr);

  {
    std::ofstream f(output_path);
    f << "Area (S),Average (S),Variance(S),Area (N),Average (N),Variance (N)\n";
    mln_foreach(auto x, tree.nodes())
      f << accu::extractor::count(res[x]) << ","
        << accu::extractor::mean(res[x]) << ","
        << accu::extractor::variance(res[x]) << ","
        << accu::extractor::count(res2[x]) << ","
        << accu::extractor::mean(res2[x]) << ","
        << accu::extractor::variance(res2[x])
        << "\n";
  }
}
