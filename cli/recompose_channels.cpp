#include <mln/core/image/image2d.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/rgba.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
  if (argc < 5)
    {
      std::cerr << "Usage: " << argv[0] << " out.tiff channel-0 channel-1...\n"
        "Recompose a color image from its individual channels...\n";
      std::exit(1);
    }

  using namespace mln;

  if (argc == 5) {
    image2d<rgb16>  out;
    image2d<uint16> f;

    for (int i = 2; i < 5; ++i) {
      io::imread(argv[i], f);
      if (i == 2)
        resize(out, f);
      copy(f, channel(out, i-2));
    }
    io::imsave(out, argv[1]);
  }
  else if (argc == 6) {
    image2d<colors::rgba16>  out;
    image2d<uint16> f;

    for (int i = 2; i < 6; ++i) {
      io::imread(argv[i], f);
      if (i == 2)
        resize(out, f);
      copy(f, channel(out, i-2));
    }
    io::imsave(out, argv[1]);
  } else {
    std::cerr << "To many channels.\n";
    std::exit(1);
  }

}
