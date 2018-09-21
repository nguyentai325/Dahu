#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/io/imread.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <apps/g2/compute_ctos.hpp>



int main(int argc, char** argv)
{
  using namespace mln;


  typedef morpho::component_tree<unsigned, image2d<unsigned> > tree_t;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " image1 [image2...]\n"
              << "Ouputs some stats about the (Color) ToS (#nodes,Average Depth,Max Depth)\n";
    std::exit(1);
  }

  image2d<rgb8> f;
  tree_t tree;

  std::cout << "Image,Size,Count,AvgDepth,MaxDepth\n";

  for (int i = 1; i < argc; ++i)
    {
      try {
        io::imread(argv[i], f);
      } catch (...) {
        image2d<uint8> g;
        io::imread(argv[i], g);
        f = eval(imcast<rgb8>(g));
      }

      tree = compute_ctos(f);

      auto depthmap = morpho::compute_depth(tree);

      size_t n = 0;
      double sumdepth = 0;
      unsigned maxdepth = 0;
      mln_foreach (auto x, tree.nodes())
        {
          sumdepth += depthmap[x];
          maxdepth = std::max(maxdepth, depthmap[x]);
          n++;
        }

      std::cout << argv[i] << "," << f.domain().size() << ","
                << n << "," << (sumdepth / n) << "," << maxdepth << "\n";
    }
}


