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
        "Convert a 8-bit Lab* image ot a 8-bit RGB image (inverse of rgb2lab)\n";
      std::exit(1);
    }

  using namespace mln;

  image2d<rgb8> f;
  io::imread(argv[1], f);


  auto h = transform(f, [](rgb8 x) -> rgb8 {
      lab<float> v = {
        x[0] * 100.0f / 256.0f,
        x[1] * 220.0f / 256.0f - 110.0f,
        x[2] * 220.0f / 256.0f - 110.0f
      };

      return lab2rgb<uint8>(v);
    });

  io::imsave(h, argv[2]);
}
