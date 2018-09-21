#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/morpho/filtering.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>
#include "topology.hpp"
#include <libgen.h>
#include "addborder.hpp"


void usage(int argc, char** argv)
{
  if (argc < 3)
    {
      std::cerr << "Usage: " << argv[0] << " ima.(pgm) out_without_extension lambda1 [lambda2  [lambda3 [...]]]" << std::endl
		<< "Compute the ToS marginally, compute area attribute and compute another ToS on it."
		<< "Output the results of filtering" << std::endl;
      abort();
    }
}




int main(int argc, char** argv)
{
  using namespace mln;

  usage(argc, argv);

  const char* filename = argv[1];

  image2d<uint8> ima;
  io::imread(filename, ima);

  typedef UInt<9> V;
  typedef image2d<V> I;
  I x2 = transform(ima, [](uint8 v) -> V { return v * 2; });
  I ori = addborder(x2);

  image2d<V> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;
  std::tie(K, parent, S) = morpho::ToS(ori, c4);

  point2d strides = point2d{2,2};

  for (int i = 3; i < argc; ++i) {
    int fvalue = std::atoi(argv[i]);
    image2d<V> out;
    image2d<uint8> out2;
    out = morpho::area_filter(K, K, parent, S, fvalue, K1::is_face_2, unsigned ());
    out2 = transform(out, [](V x) -> uint8 { return x / 2; });

    image2d<uint8> under;
    resize(under, ori);
    copy(out2 | sbox2d(out2.domain().pmin, out2.domain().pmax, strides), under);
    io::imsave(out2, boost::str(boost::format("%s-%06i_k1.tiff") % argv[2] % fvalue).c_str());
    io::imsave(under, boost::str(boost::format("%s-%06i.tiff") % argv[2] % fvalue).c_str());
  }

}
