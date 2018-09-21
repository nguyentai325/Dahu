#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include "curvature.hpp"


int main(int argc, char** argv)
{
  if (argc != 3)
    {
      std::cerr << "Usage: " << argv[0] << "input.ppm output.tiff" << std::endl; 
      std::terminate();
    }

  using namespace mln;

  image2d<rgb8> input;
  io::imread(argv[1], input);

  image2d<float> out = curvature(input);
  io::imsave(out, argv[2]);

}
