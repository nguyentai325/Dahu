#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/copy.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>
#include <stack>
#include "topology.hpp"
#include "addborder.hpp"

namespace mln
{

  template <typename V>
  void
  filter_setmask(const image2d<V>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S, unsigned grain, image2d<uint8>& mask, uint8 lbl)
  {
    image2d<unsigned> area;
    resize(area, K);
    for (int i = S.size() - 1; i >= 0; --i)
      {
	unsigned p = S[i];
	unsigned rp = K[p] == K[parent[p]] ? parent[p] : p;

	if (K1::is_face_2(K.point_at_index(p)))
	  ++area[rp];
      }

    for (unsigned p: S)
      {
	unsigned rp = K[p] == K[parent[p]] ? parent[p] : p;
	if (area[rp] < grain)
	  mask[p] |= lbl;
      }
  }



  image2d<rgb8>
  setmask_with_mean(const image2d<rgb8>& ima, const image2d<uint8>& mask)
  {
    static constexpr unsigned PROCESSED = value_traits<unsigned>::max();
    static constexpr unsigned UNPROCESSED = 0;

    image2d<unsigned> lbl;
    resize(lbl, mask, mask.border(), 0);
    extension::fill(lbl, PROCESSED);

    unsigned sz = mask.domain().size();

    std::vector<unsigned> queue_;
    queue_.reserve(sz);

    std::stack<unsigned, std::vector<unsigned> > q(std::move(queue_));
    std::vector< rgb<unsigned> > sum_(sz+1, rgb<unsigned>{0,0,0} );
    std::vector<unsigned>	 count_(sz+1, 0);


    mln_pixter(px, lbl);
    auto didx = wrt_delta_index(mask, c4.dpoints);

    unsigned nlabel = 1;

    {
      mln_forall(px)
	if (px->val() == UNPROCESSED)
	  {
	    unsigned i = px->index();
	    uint8 v = mask[i];
	    if (mask[i] == 0) {
	      lbl[i] = 0;
	      continue;
	    }
	    q.push(i);
	    while (not q.empty())
	      {
	    	unsigned p = q.top();
	    	q.pop();
	    	lbl[p] = nlabel;

	    	point2d pp = mask.point_at_index(p);
	    	if (K1::is_face_2(pp))
	    	  {
	    	    sum_[nlabel] += ima(pp/2);
	    	    count_[nlabel]++;
	    	  }
	    	mln_foreach (int k, didx)
	    	  {
	    	    unsigned n = p+k;
	    	    if (lbl[n] == UNPROCESSED and mask[n] == v)
	    	      q.push(n);
	    	  }
	    }
	    ++nlabel;
	  }
    }
    std::cout << "Nlabel: " << nlabel << std::endl;
    assert(nlabel < (sz+1));

    image2d<rgb8> out;
    resize(out, mask, mask.border(), rgb8{0,0,255});
    mln_pixter(pin, pout, lbl, out);
    mln_forall(pin, pout)
      {
	if (pin->val() == 0 and K1::is_face_2(pin->point()))
	  pout->val() = ima(pin->point() / 2);
	else if (count_[pin->val()] > 0)
	  pout->val() = sum_[pin->val()] / count_[pin->val()];
      }

    return out;

  }
}

void usage(int argc, char** argv)
{
  if (argc < 4)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm out_wo_ext grain1 [grain2 [grain3...]]" << std::endl;
      std::abort();
    }
}



int main(int argc, char** argv)
{
  using namespace mln;

  usage(argc, argv);

  std::string filename = argv[1];
  std::string stem = argv[2];

  image2d<rgb8> ima;
  io::imread(filename, ima);

  typedef UInt<9> V;
  typedef image2d<V> I;
  I r = transform(ima, [](rgb8 v) -> V { return v[0] * 2; });
  I g = transform(ima, [](rgb8 v) -> V { return v[1] * 2; });
  I b = transform(ima, [](rgb8 v) -> V { return v[2] * 2; });
  I rr = addborder(r);
  I gg = addborder(g);
  I bb = addborder(b);


  image2d<V> rK, gK, bK;
  image2d<unsigned> rparent, gparent, bparent;
  std::vector<unsigned> rS, gS, bS;

  std::tie(rK, rparent, rS) = morpho::ToS(rr, c4);
  std::tie(gK, gparent, gS) = morpho::ToS(gg, c4);
  std::tie(bK, bparent, bS) = morpho::ToS(bb, c4);

  auto bima = addborder(ima);
  point2d strides{2,2};
  for (int i = 3; i < argc; ++i)
    {
      unsigned grain = std::atoi(argv[i]);

      image2d<uint8> mask;
      resize(mask, rK, rK.border(), 0);
      filter_setmask(rK, rparent, rS, grain, mask, 0x04);
      filter_setmask(gK, gparent, gS, grain, mask, 0x02);
      filter_setmask(bK, bparent, bS, grain, mask, 0x01);
      io::imsave(mask, boost::str(boost::format("%s-%06i.tiff") % "mask" % grain).c_str());

      auto out = setmask_with_mean(bima, mask);
      io::imsave(out, boost::str(boost::format("%s-%06i_k1.tiff") % stem % grain).c_str());

      image2d<rgb8> under;
      resize(under, bima);
      copy(out | sbox2d(out.domain().pmin, out.domain().pmax, strides), under);
      io::imsave(under, boost::str(boost::format("%s-%06i.tiff") % stem % grain).c_str());
    }
}
