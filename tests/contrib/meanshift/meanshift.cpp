#include <mln/contrib/meanshift/meanshift.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

int main(int argc, const char** argv)
{
  if (argc < 5)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm hs hr output.hpp\n";
      std::exit(1);
    }

  using namespace mln;

  image2d<rgb8> f;

  io::imread(argv[1], f);

  float hs = std::atof(argv[2]);
  float hr = std::atof(argv[3]);
  image2d<rgb8> out = contrib::meanshift(f, hs, hr);

  io::imsave(out, argv[4]);
}
