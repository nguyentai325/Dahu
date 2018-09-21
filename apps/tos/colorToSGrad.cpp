#include "colorToSGrad.hpp"

#include <mln/core/algorithm/transform.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/morpho/maxtree_pqueue_parallel.hpp>
#include <mln/morpho/filtering.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/gradient.hpp>
#include <apps/tos/topology.hpp>
#include <apps/tos/Kinterpolate.hpp>

namespace mln
{

  void colorToSGrad(const image2d<rgb8>& ima,
		    image2d<unsigned>& K,
		    image2d<unsigned>& parent,
		    std::vector<unsigned>& S)
  {
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

    auto r_area = morpho::area_compute(rK, rparent, rS, K1::is_face_2);
    auto g_area = morpho::area_compute(gK, gparent, gS, K1::is_face_2);
    auto b_area = morpho::area_compute(bK, bparent, bS, K1::is_face_2);

    int size = 7;
    auto grad_r = interpolate_k1(gradient(rr, size));
    auto grad_g = interpolate_k1(gradient(gg, size));
    auto grad_b = interpolate_k1(gradient(bb, size));

    image2d<unsigned> area;
    resize(area, r_area);
    mln_foreach(const point2d& p, area.domain())
      {
	int r = grad_r(p), g = grad_g(p), b = grad_b(p);
	int ar = r_area(p), ag = g_area(p), ab = b_area(p);

	if (r >= g) {
	  area(p) = (r >= b) ? ar : ab;
	} else {
	  area(p) = (g >= b) ? ag : ab;
	}
      }

    std::tie(K, parent, S) = morpho::ToS_priority(area, c4);
  }

  void colorToSGrad_2f(const image2d<rgb8>& ima,
		       image2d<unsigned>& K,
		       image2d<unsigned>& parent,
		       std::vector<unsigned>& S)
  {
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

    auto r_area = morpho::area_compute(rK, rparent, rS, K1::is_face_2);
    auto g_area = morpho::area_compute(gK, gparent, gS, K1::is_face_2);
    auto b_area = morpho::area_compute(bK, bparent, bS, K1::is_face_2);

    int size = 7;
    auto grad_r = gradient(rr, size);
    auto grad_g = gradient(gg, size);
    auto grad_b = gradient(bb, size);

    image2d<unsigned> area;
    resize(area, rr);
    mln_foreach(const point2d& p, area.domain())
      {
	int r = grad_r(p), g = grad_g(p), b = grad_b(p);
	int ar = r_area(2*p), ag = g_area(2*p), ab = b_area(2*p);

	if (r >= g) {
	  area(p) = (r >= b) ? ar : ab;
	} else {
	  area(p) = (g >= b) ? ag : ab;
	}
      }

    std::tie(K, parent, S) = morpho::ToS(area, c4);
  }


 void colorToSGrad_with_mintree_(const image2d<rgb8>& ima,
				 image2d<unsigned>& K,
				 image2d<unsigned>& parent,
				 std::vector<unsigned>& S)
  {
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

    auto r_area = morpho::area_compute(rK, rparent, rS, K1::is_face_2);
    auto g_area = morpho::area_compute(gK, gparent, gS, K1::is_face_2);
    auto b_area = morpho::area_compute(bK, bparent, bS, K1::is_face_2);

    int size = 7;
    auto grad_r = interpolate_k1(gradient(rr, size));
    auto grad_g = interpolate_k1(gradient(gg, size));
    auto grad_b = interpolate_k1(gradient(bb, size));

    image2d<unsigned> area;
    resize(area, r_area);
    mln_foreach(const point2d& p, area.domain())
      {
	int r = grad_r(p), g = grad_g(p), b = grad_b(p);
	int ar = r_area(p), ag = g_area(p), ab = b_area(p);

	if (r >= g) {
	  area(p) = (r >= b) ? ar : ab;
	} else {
	  area(p) = (g >= b) ? ag : ab;
	}
      }

    std::tie(parent, S) = morpho::impl::serial::maxtree_pqueue(area, c4, std::greater<unsigned> ());
    K = area;
  }

 void colorToSGrad_with_mintree(const image2d<rgb8>&	ima,
				image2d<unsigned>&	K,
				image2d<unsigned>&	parent,
				std::vector<unsigned>&	S)
  {
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

    auto r_area = morpho::area_compute(rK, rparent, rS, K1::is_face_2);
    auto g_area = morpho::area_compute(gK, gparent, gS, K1::is_face_2);
    auto b_area = morpho::area_compute(bK, bparent, bS, K1::is_face_2);

    int size = 7;
    auto grad_r = interpolate_k1(gradient(rr, size));
    auto grad_g = interpolate_k1(gradient(gg, size));
    auto grad_b = interpolate_k1(gradient(bb, size));

    image2d<unsigned> area;
    resize(area, r_area).init(0);

    point2d p;
    mln_iter(q, c8(p));
    sbox2d dom_2f = { area.domain().pmin, area.domain().pmax, point2d{2,2} };
    mln_foreach(p, dom_2f)
      {
	int r = grad_r(p), g = grad_g(p), b = grad_b(p);
	int ar = r_area(p), ag = g_area(p), ab = b_area(p);

	if (r >= g) {
	  area(p) = (r >= b) ? ar : ab;
	} else {
	  area(p) = (g >= b) ? ag : ab;
	}

	mln_forall(q)
	  area.at(*q) = std::max(area.at(*q), area(p));
      }

    std::tie(parent, S) = morpho::impl::serial::maxtree_pqueue(area, c4, std::greater<unsigned> ());
    K = area;
  }


} // end of namespace mln
