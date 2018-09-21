#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/colors/lsh.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/labeling/blobs.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>

#include <mln/graphcut/graphcut.hh>

#include <mln/morpho/tos/tos2.hpp>
#include "addborder.hpp"
#include "thicken.hpp"
#include "topology.hpp"

namespace mln
{

  template <typename V, typename T, class FilterFun>
  image2d<V>
  setmean_on_nodes(const image2d<V>& ima, const image2d<T>& K,
		   const image2d<unsigned>& parent, const std::vector<unsigned>& S, FilterFun pred)
  {
    typedef rgb<unsigned> SumType;
    image2d<SumType>  sum;
    image2d<unsigned> count;
    resize(count, K).init(0);
    resize(sum, K).init(SumType ());

    for (int i = S.size() - 1; i >= 0; --i)
      {
	unsigned p = S[i];
	if (pred(K.point_at_index(p)))
	  {
	    unsigned q = parent[p];
	    if (K[p] != K[q])
	      q = p;
	    count[q] += 1;
	    sum[q] += ima[p];
	  }
      }

    image2d<V> out;
    resize(out, ima);
    unsigned p = S[0];
    out[p] = sum[p] / count[p];
    for (unsigned p: S)
      {
	unsigned q = parent[p];
	if (K[p] == K[q])
	  out[p] = out[q];
	else
	  {
	    while (count[p] == 0)
	      p = parent[p];
	    out[p] = sum[p] / count[p];
	  }
      }
    return out;
  }

  template <typename V>
  image2d<V>
  k1tok0(const image2d<V>& ima)
  {
    image2d<V> out(ima.nrows() / 2 + 1, ima.ncols() / 2 + 1);
    auto k1tok0dom = sbox2d(ima.domain().pmin, ima.domain().pmax, point2d{2,2});
    copy(ima | k1tok0dom, out);
    return out;
  }

}


int main(int argc, const char** argv)
{
  using namespace mln;
  typedef UInt<9> uint9;

  if (argc < 3)
    {
      std::cout << "Usage: " << argv[0] << " input output(wo ext) [saturation]" << std::endl;
      abort();
    }

  std::string filename = argv[1];
  std::string output = argv[2];
  int sat = (argc > 3) ? std::atoi(argv[3]) : 45;

  image2d<rgb8> ima;
  io::imread(filename, ima);

  image2d<lsh8> f = transform(ima, rgb2lsh<uint8> );


  io::imsave(transform(f, [sat] (lsh8 x) -> bool { return x[1] > sat; }), (boost::format("%s-mask.tiff") % output).str().c_str());
  io::imsave(transform(f, [sat] (lsh8 x) -> uint8 { return x[1]; }), (boost::format("%s-sat.tiff") % output).str().c_str());
  io::imsave(transform(f, [sat] (lsh8 x) -> uint8 { return x[1] < sat ? x[0] : 0; }), (boost::format("%s-lum.tiff") % output).str().c_str());
  io::imsave(transform(f, [sat] (lsh8 x) -> uint8 { return x[1] > sat ? x[2] : 0; }), (boost::format("%s-hue.tiff") % output).str().c_str());

  image2d<rgb8> bima = addborder(ima);
  image2d<rgb8> tmp(2*bima.nrows()-1, 2*bima.ncols()-1);
  fill(tmp, rgb8{0,0,255});
  copy(bima, tmp | sbox2d(tmp.domain().pmin, tmp.domain().pmax, point2d{2,2}));

  image2d<rgb8> final;
  resize(final, bima);
  fill(final, rgb8{0,0,255});

  //image2d<bool> mask_ = transform(f, [sat] (lsh8 x) { return (x[1] < sat); });
  image2d<bool> mask_;
  resize(mask_, f);

  // regularisation
  graphcut::graphcut(f, mask_, c4,
		     [sat] (bool x, const lsh8& v) -> float { return 0 * (x != (v[1] > sat)) + 1 * std::abs(x - 1.0 / (1.0 + std::exp(-(v[1] - sat)/3.0))); },
		     [] (bool x, bool y) { return 1 * (x != y); },
  		     [] (float x, float y) { return (std::abs(x - y) > 0.001) * (x - y); }  // to compare float
		     );


  io::imsave(mask_, (boost::format("%s-mask-reg.tiff") % output).str().c_str());

  //std::exit(1);

  // Compute TOS on cc[L]
  {
    //image2d<bool> mask = transform(f, [sat] (lsh8 x) { return (x[1] < sat); });
    image2d<bool> mask = mask_;
    auto L = transform(f, [] (lsh8 x) -> uint9 { return 2 * x[0]; });

    image2d<uint16> lbl;
    unsigned nlabel;
    //std::cout << "Start: labeling" << std::endl;
    std::tie(lbl, nlabel) = labeling::blobs(mask, c4, uint16());
    //std::cout << "End: labeling" << std::endl;
    //std::vector< std::pair< FIXME, FIXME >

    for (unsigned i = 1; i <= nlabel; ++i)
      {
	//auto x = (lbl == i);
	//auto subima = (L | ));
	image2d<uint9> ima2;
	image2d<bool> mask2;

	std::tie(ima2, mask2) = addborder2(L, lbl == i);
	//io::imsave(transform(ima2, [] (uint9 x) -> uint8 { return x/2; }),
	//	   (boost::format("L_%02i.tiff") % i).str().c_str());
	auto domain = rng::filter(ima2.domain(), [mask2](point2d p) { return mask2(p); });

	image2d<uint9> K;
	image2d<unsigned> parent;
	std::vector<unsigned> S;

	//std::cout << "Start: subtos (" << i << ")" << std::endl;
	std::tie(K, parent, S) = morpho::subToS(ima2, mask2, domain, c4);
	//std::cout << "End: subtos (" << i << ")" << std::endl;

	//std::cout << "Start: thicken (" << i << ")" << std::endl;
	auto w = thicken_tdn(K, parent, S, c8);
	//std::cout << "End: thicken (" << i << ")" << std::endl;

	auto res = setmean_on_nodes(tmp, w, parent, S, K1::is_face_2);
	auto k0 = k1tok0(res);
	//io::imsave(res, (boost::format("Lres_%02i.tiff") % i).str().c_str());
	copy(k0 | mask2, final | mask2);
	//fill(final|mask2, rgb8{0,255,0});
      }
  }

  // Compute TOS on cc[H]
  {
    auto mask = lnot(mask_);
    //image2d<bool> mask = transform(f, [sat] (lsh8 x) { return (x[1] >= sat); });
    auto H = transform(f, [] (lsh8 x) -> uint9 { return 2 * x[2]; });

    image2d<uint16> lbl;
    unsigned nlabel;
    std::tie(lbl, nlabel) = labeling::blobs(mask, c8, uint16());

    for (unsigned i = 1; i <= nlabel; ++i)
      {
	image2d<uint9> ima2;
	image2d<bool> mask2;

	std::tie(ima2, mask2) = addborder2(H, lbl == i);
	//io::imsave(transform(ima2, [] (uint9 x) -> uint8 { return x/2; }),
	//	   (boost::format("H_%02i.tiff") % i).str().c_str());
	auto domain = rng::filter(ima2.domain(), [mask2](point2d p) { return mask2(p); });

	image2d<uint9> K;
	image2d<unsigned> parent;
	std::vector<unsigned> S;

	std::tie(K, parent, S) = morpho::subToS(ima2, mask2, domain, c4);
	auto w = thicken_tdn(K, parent, S, c8);

	auto res = setmean_on_nodes(tmp, w, parent, S, K1::is_face_2);
	//io::imsave(res, (boost::format("Hres_%02i.tiff") % i).str().c_str());
	copy(k1tok0(res) | mask2, final | mask2);
	//fill(final|mask2, rgb8{0,0,255});
      }
  }

  io::imsave(final, (boost::format("%s.tiff") % output).str().c_str() );
}
