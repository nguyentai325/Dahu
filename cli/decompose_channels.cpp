#include <mln/core/image/image2d.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/rgba.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
  if (argc != 3)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm out.tiff\n"
        "Decompose a color image intro individual channel out-0.tiff out-1.tiff out-2.tiff...\n";
      std::exit(1);
    }

  using namespace mln;

  std::string oname = argv[2];
  int pos = oname.rfind('.');
  std::string fmt = oname.substr(0, pos) + "-%i" + oname.substr(pos);

  try {
    image2d<colors::rgba8> f;
    io::imread(argv[1], f);

    for (int i = 0; i < 4; ++i)
      io::imsave(channel(f,i), (boost::format(fmt) % i).str());
  } catch (...) {
    image2d<rgb8> f;
    io::imread(argv[1], f);

    for (int i = 0; i < 3; ++i)
      io::imsave(channel(f,i), (boost::format(fmt) % i).str());
  }

}
