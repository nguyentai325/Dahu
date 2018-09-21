#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/algorithm/transform.hpp>
#include <mln/morpho/algebraic_filter.hpp>
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

  image2d<uint8> out;
  unsigned lambda = std::atoi(argv[3]);

  if (argv[2][0] == '4')
    out = morpho::area_closing(ima, c4, lambda);
  else
    out = morpho::area_closing(ima, c8, lambda);

  io::imsave(out, argv[4]);
}
