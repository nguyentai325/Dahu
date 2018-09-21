#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/algorithm/transform.hpp>
#include "gradient.hpp"

void usage(char* argv[])
{
  std::cerr << "usage: " << argv[0] << " input.ppm size output.ppm" << std::endl;
  std::abort();
}


int main(int argc, char** argv)
{
  if (argc != 4)
    usage(argv);

  using namespace mln;


  image2d<rgb8> ima;
  io::imread(argv[1], ima);

  int size = std::atoi(argv[2]);
  if (size < 3 || (size % 2 == 0))
    usage(argv);

  typedef uint8 V;
  auto r = transform(ima, [](rgb8 v) -> V { return v[0]; });
  auto g = transform(ima, [](rgb8 v) -> V { return v[1]; });
  auto b = transform(ima, [](rgb8 v) -> V { return v[2]; });

  auto gr = gradient(r, size);
  auto gg = gradient(g, size);
  auto gb = gradient(b, size);

  io::imsave(gr, "gr.tiff");
  io::imsave(gg, "gg.tiff");
  io::imsave(gb, "gb.tiff");


  image2d<rgb8> output = transform(imzip(gr, gg, gb), [] (const std::tuple<V,V,V>& x) {
      V r,g,b;
      std::tie(r,g,b) = x;
      if (r >= g) {
	return (r >= b) ? rgb8{255,0,0} : rgb8{0,0,255};
      } else {
	return (g >= b) ? rgb8{0,255,0} : rgb8{0,0,255};
      }
    });


  io::imsave(output, argv[3]);

}




