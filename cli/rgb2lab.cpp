#include <mln/core/image/image2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/colors.hpp>
#include <mln/colors/lab.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

int main(int argc, char** argv)
{
  if (argc != 3)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm out.tiff\n"
        "Convert a 8-bit rgb image ot a 8-bit La*b* image\n"
        "L is in [0-100] and a*b* in [0-110], the values are stetched"
        "linearly in the range [0-255]";
      std::exit(1);
    }

  using namespace mln;

  image2d<rgb8> f;
  io::imread(argv[1], f);


  auto h = transform(f, [](rgb8 x) -> rgb8 {
      lab<float> y = rgb2lab(x);
      return rgb8{ uint8(y[0] * 256 / 100),
          uint8( (y[1] + 110) * 256 / 220),
          uint8( (y[2] + 110) * 256 / 220) };
    });


  io::imsave(h, argv[2]);
}
