#include <mln/core/image/image2d.hpp>
#include <mln/io/imread.hpp>
#include <apps/tos/addborder.hpp>

#include "dispsaliency.hpp"
#include <mln/io/imsave.hpp>

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage:" << argv[0] << " saliency.tiff [original.tiff]" << std::endl
	      << "Saliency must be a saliency map (n x m) with saliency on edges." << std::endl
	      << "Original must be either an rgb or gray level image with size (n-3)/2 x (m-3)/2 (wo border)"
	      << std::endl;
    std::terminate();
  }


  QApplication a(argc, argv);
  const char* s_input = argv[1];
  const char* s_input_2 = NULL;

  if (argc == 3)
    s_input_2 = argv[2];


  image2d<float> smap;
  io::imread(s_input, smap);



  if (s_input_2 != NULL)
    {
      image2d<rgb8> ref;
      try {
	io::imread(s_input_2, ref);
      } catch (io::MLNIOException) {
	image2d<uint8> ref_;
	io::imread(s_input_2, ref_);
	ref = transform(ref_, [](uint8 x) { return rgb8{x,x,x}; });
      }
      ref = addborder(ref, lexicographicalorder_less<rgb8>() );

      point2d sz = smap.domain().shape();
      point2d sz2 = ref.domain().shape();
      if (sz != (sz2 * 2 - 1))
	throw std::runtime_error("Domains have invalid size");

      Displayer disp(smap, &ref, argv[1]);
      a.exec();
    }
  else
    {
      Displayer disp(smap, NULL, argv[1]);
      a.exec();
    }

}
