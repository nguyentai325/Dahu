# include <mln/core/image/image2d.hpp>

# include <mln/io/imread.hpp>
# include <mln/io/imsave.hpp>

# include <mln/core/algorithm/paste.hpp>
# include <mln/core/algorithm/transform.hpp>


#include <mln/morpho/tos/tos.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/mumford_shah.hpp>
#include <apps/tos/addborder.hpp>
#include <apps/tos/set_mean_on_nodes.hpp>
#include <apps/tos/objdetection.hpp>
#include <boost/format.hpp>



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input(gray) output[wo ext] [lambda [lambda_2...]]" << std::endl
	    << "Perform a simplification of the ToS by removing non significant"
	    << "level lines. A closing in shape space in performed." << std::endl
	    << "Param:" << std::endl
	    << "lambda: minimal grain size (in shape space) (default: 0)" << std::endl
	    << "lambda_i: ..." << std::endl
	    << "If a single lambda is given, the output image is named <output>.tiff, otherwise "
	    << "it creates as many images as lambdas named <output>-<lambda_i>.tiff"
	    << std::endl;
  std::terminate();
}

namespace mln
{

  template <typename V>
  image2d<unsigned>
  make_unique_id(const image2d<V>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    image2d<unsigned> out;
    resize(out, K);
    for (unsigned x: S)
      if (K[x] == K[parent[x]])
	out[x] = parent[x];
      else
	out[x] = x;
    return out;
  }

}


int main(int argc, char** argv)
{
  if (argc < 4)
    usage(argv);

  using namespace mln;
  image2d<uint8> ima;
  io::imread(argv[1], ima);

  typedef UInt<9> V;
  image2d<uint8> bima = addborder(ima);
  image2d<V> ima_ = transform(bima, [](uint8 x) -> V { return x*2; });
  image2d<uint8> f = interpolate_k1(bima);


  image2d<V> K;
  std::vector<unsigned> S;
  image2d<unsigned> parent;

  std::tie(K, parent, S) = morpho::ToS(ima_, c4);
  image2d<float> energy = compute_energy(bima, K, parent, S);
  energy[S[0]] = value_traits<float>::max();

  image2d<unsigned> ids = make_unique_id(K, parent, S);

  image2d<uint8> out, final;
  resize(final, ima);
  resize(out, K);
  for (int i = 3; i < argc; ++i)
    {
      int lambda = std::atoi(argv[i]);
      auto cIds = clone(ids);
      auto cParent = clone(parent);
      auto cEnergy = clone(energy);

      getObjects(cEnergy, cIds, cParent, S, lambda);
      out = set_mean_on_node(f, cIds, S, cParent, K1::is_face_2);

      // // recontruct
      // {
      // 	for (unsigned x: S)
      // 	  out[x] = (float) cIds[x] / out.domain().size();
      // }

      // Subsample
      {
	box2d d = out.domain();
	sbox2d sub_domain { {2,2}, {d.pmax[0]-2, d.pmax[1]-2}, {2,2} };
	copy(out | sub_domain, final);
      }

      if (argc == 4)
	io::imsave(final, (std::string(argv[2]) + ".tiff").c_str());
      else
	io::imsave(final, (boost::format("%s-%06i.tiff") % argv[2] % lambda).str().c_str());
    }
}
