# include <mln/core/image/image2d.hpp>

# include <mln/io/imread.hpp>
# include <mln/io/imsave.hpp>

# include <mln/core/algorithm/copy.hpp>
# include <mln/core/algorithm/transform.hpp>
# include <mln/core/algorithm/accumulate.hpp>
# include <mln/accu/accumulators/minmax.hpp>


#include <mln/morpho/tos/tos.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/addborder.hpp>
#include <apps/tos/routines.hpp>
//#include <apps/tos/objdetection.hpp>

#include <apps/attributes/meaningfullness.hpp>
#include <apps/saliency/extinction.hpp>
#include <apps/saliency/saliency.hpp>
#include <apps/llview/llview.hpp>

#include <boost/format.hpp>



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input(gray) output[wo ext] α β γ ε grain [λ₁ [λ₂...]]" << std::endl
	    << "Perform a simplification of the ToS by removing non significant"
	    << "level lines. A closing in shape space in performed." << std::endl
	    << "Param:" << std::endl
	    << "  α, β, γ, ε:	meaningfullness energy parameters" << std::endl
	    << "  grain:	grain filter (in image domain)" << std::endl
	    << "  λ₁, λ₂:	minimal dynamic (in shape space) (default: 0)" << std::endl
	    << "It creates as many images as lambdas named <output>-<λ>.tiff"
	    << "The meaninfullness energy is: α.A + B + β.C with" << std::endl
	    << "   A = Context Energy Distance f(Internal Region, External Region)" << std::endl
	    << "   B = Curvature" << std::endl
	    << "   C = Constraint" << std::endl
	    << "Suggested value: "
	    << "   α: 2, β: 1, ε: 5, γ: 0.05"
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
  if (argc < 9)
    usage(argv);

  const char* infname = argv[1];
  const char* outfname = argv[2];
  float alpha = std::atof(argv[3]);
  float beta = std::atof(argv[4]);
  float gamma = std::atof(argv[5]);
  int eps = std::atof(argv[6]);
  int grain = std::atoi(argv[7]);

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

  if (grain != 0) {
    grainfilter_inplace(K, parent, S, grain, K1::is_face_2);
  }


  image2d<float> energy = meaningfullness(bima, K, parent, S,
					  alpha, beta, gamma, eps);

  //energy[S[0]] = value_traits<float>::max();
  //image2d<unsigned> ids = make_unique_id(K, parent, S);


  // Retrieve the extinction values for each node.
  // Non minima nodes are set to 0
  auto extmap = extinction(energy, K, parent, S);
  auto sal = saliencymap(extmap, K, parent, S);

  io::imsave(sal, (boost::format("%s-saliency.tiff") % argv[2]).str().c_str());


  auto realnodes = get_real_nodes(K, parent, S);
  image2d<bool> mask;
  resize(mask, K).init(false);
  for (unsigned x : realnodes)
    mask[x] = true;

  std::pair<float, float> en_minmax  = accumulate(energy | mask, accu::features::minmax<> ());
  std::pair<float, float> ext_minmax = accumulate(extmap | mask, accu::features::minmax<> ());

  std::cout << "Some statistics:" << std::endl
	    << "Energy: " << en_minmax.first  << "," << en_minmax.second << std::endl
	    << "Extinction: " << ext_minmax.first << "," << ext_minmax.second << std::endl;


  image2d<uint8> final;
  resize(final, ima);
  for (int i = 8; i < argc; ++i)
    {
      float lambda = std::atof(argv[i]);

      auto cK = clone(K);
      auto cParent = clone(parent);

      // Filter out the minima that are not enough meaningfull
      attributefilter_inplace(extmap, cK, cParent, S, lambda);


      image2d<int> out = compute_mean_on_node(f, cK, cParent, S, K1::is_face_2);

      // // recontruct
      // {
      // 	for (unsigned x: S)
      // 	  out[x] = (float) cIds[x] / out.domain().size();
      // }

      // Subsample
      {
	box2d d = out.domain();
	sbox2d sub_domain { d.pmin + 2 , d.pmax - 2, {2,2} };
	copy(out | sub_domain, final);
      }

      image2d<bool> llines = llview( eval(extmap >= lambda), K, parent, S);

      io::imsave(final, (boost::format("%s-simp-%04f.tiff") % argv[2] % lambda).str().c_str());
      io::imsave(llines, (boost::format("%s-llines-%04f.tiff") % argv[2] % lambda).str().c_str());
    }
}
