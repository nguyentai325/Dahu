#include <mln/core/image/image2d.hpp>
#include <mln/core/image/sub_image.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

void usage(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr << "Usage: " << argv[0] << "in.tiff out.tiff" << std::endl
		<< "K1 to K2 conversion." << std::endl;
      abort();
    }
}


int main(int argc, char** argv)
{
  using namespace mln;
  usage(argc, argv);


  image2d<rgb8> ima, k0;
  io::imread(argv[1], ima);

  box2d d1 = ima.domain();
  sbox2d d0(d1.pmin, d1.pmax, point2d{2,2});
  box2d d = {{0,0}, d0.shape()};
  k0.resize(d);
  copy(ima | d0, k0);

  io::imsave(k0, argv[2]);
}
