#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/morpho/closing_by_reconstruction.hpp>

#include <iostream>

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 5)
    {
      std::cerr << "Usage: " << argv[0] << " input(gray) markers(gray) connectivity output(gray)"
                << std::endl
                << "connectivity: 4 / 8" << std::endl;
      std::exit(1);
    }

  image2d<uint8> ima;
  image2d<uint8> markers;
  image2d<uint8> out;
  io::imread(argv[1], ima);
  io::imread(argv[2], markers);

  if (argv[3][0] == '4')
    out = morpho::closing_by_reconstruction(ima, markers, c4);
  else
    out = morpho::closing_by_reconstruction(ima, markers, c8);

  io::imsave(out, argv[4]);
}
