#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/morpho/structural/dilate.hpp>
#include <mln/core/win2d.hpp>
#include <mln/core/se/ball.hpp>
#include <iostream>

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 4)
    {
      std::cerr << "Usage: " << argv[0] << " input ball_radius output(gray)" << std::endl;
      std::exit(1);
    }

  image2d<uint8> ima;
  io::imread(argv[1], ima);

  int radius = std::atoi(argv[2]);
  se::ball2d w = se::ball2d::make(radius + 0.3f);

  image2d<uint8> out = morpho::structural::dilate(ima, w);
  io::imsave(out, argv[3]);
}
