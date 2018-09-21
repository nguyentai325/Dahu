#ifndef MUMFORD_SHAH_HPP
# define MUMFORD_SHAH_HPP

# include <mln/core/trace.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/core/neighb2d.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/io/imprint.hpp>
# include "topology.hpp"

namespace mln
{

  /// \brief Take an image and puts its curvature on the edges
  /// The image is doubled, non-edges pixels have undefined values.
  template <typename V>
  image2d<float>
  curvature_on_edge(const image2d<V>& ima);


  /// \brief Compute the mumford-shah energy for each nodes.
  /// \pre K must be twice as big as ima.
  template <typename V, typename T>
  image2d<float>
  compute_energy(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent,
		 const std::vector<unsigned>& S, int eps = 5);


  /// \brief Perform a closing on an attribute image
  template <typename V>
  image2d<float>
  close(const image2d<float>& energy, const image2d<V>& K, const image2d<unsigned>& parent,
	const std::vector<unsigned>& S);

  /****************************/
  /** Implementation         **/
  /****************************/

  namespace internal
  {
    template <typename V>
    decltype( std::declval<V>() * 2 )
      sqr(const V& x)
    {
      return x*x;
    }

    template <typename V>
    typename std::enable_if< std::is_scalar<V>::value, V>::type
    sum(const V& x)
    {
      return x;
    }

    template <typename V>
    auto
    sum(const V& x) -> typename std::enable_if< !std::is_scalar<V>::value, decltype(x[0] + x[0])>::type
    {
      typedef decltype( x[0] + x[0] ) R;
      R s = x[0];
      for (int i = 1; i < V::ndim; ++i)
	s += x[i];
      return s;
    }

    template <typename V, bool scalar = std::is_scalar<V>::value >
    struct energy_helper
    {
      typedef uint64 sum_type;
      typedef uint64 sum_sqr_type;
    };

    template <typename V>
    struct energy_helper<V, false>
    {
      typedef vec<uint64, V::ndim> sum_type;
      typedef vec<uint64, V::ndim> sum_sqr_type;
    };

    template <typename V>
    struct energy_t
    {
      typedef typename energy_helper<V>::sum_type sum_type;
      typedef typename energy_helper<V>::sum_sqr_type sum_sqr_type;

      static constexpr const float alpha = 1;
      static constexpr const float beta = 100;

      bool dejavu = false;
      int depth = -1;
      int m_e_length = 0;
      float m_e_sumcurv = 0;
      int m_v_n_int = 0;
      sum_type     m_v_sum_int = literal::zero;
      sum_sqr_type m_v_sum_int_sqr = literal::zero;
      int m_v_n_ext = 0;
      sum_type     m_v_sum_ext = literal::zero;
      sum_sqr_type m_v_sum_ext_sqr = literal::zero;

      energy_t()
      {
      }


      float external_energy() const {
	float Vin = sum(m_v_sum_int_sqr) - sum(sqr(m_v_sum_int)) / m_v_n_int;
	float Vout = m_v_n_ext == 0 ? 0 : sum(m_v_sum_ext_sqr) - sum(sqr(m_v_sum_ext)) / m_v_n_ext;
	float Eext;
	if (Vin + Vout == 0)
	  Eext = 0;
	else
	  Eext = (Vin + Vout) / (sum(m_v_sum_int_sqr + m_v_sum_ext_sqr) -
				 sum(sqr(m_v_sum_int + m_v_sum_ext)) / (m_v_n_ext + m_v_n_int));

	return Eext;
      }

      float energy() const
      {
	//auto sqr = [](float x) { return x*x; };
	float Eext = external_energy();
	float Eint = m_e_sumcurv / m_e_length;
	float Econ = 1.0 / m_e_length;
	// std::cout << "Contour length: " <<  m_e_length << std::endl
	// 	  << "Curvature: " << m_e_sumcurv << std::endl;
	assert(m_e_length > 0);
	assert(Eext >= 0);
	assert(Eint >= 0);
	assert(Econ >= 0);
	return alpha * Eint + Eext + beta * Econ;
      }

      friend
      std::ostream&
      operator <<(std::ostream& os, const energy_t& acc)
      {
	return os << "Contour length: " <<  acc.m_e_length << std::endl
		  << "Curvature: " << acc.m_e_sumcurv << std::endl
		  << "Internal Vertex: " << acc.m_v_n_int << std::endl
		  << "External Vertex: " << acc.m_v_n_ext << std::endl
		  << "Internal Sum / Sum**2: " << acc.m_v_sum_int << " / " << acc.m_v_sum_int_sqr << std::endl
		  << "External Sum / Sum**2: " << acc.m_v_sum_ext << " / " << acc.m_v_sum_ext_sqr << std::endl
		  << "Internal Variance: " << (sum(acc.m_v_sum_int_sqr) - sum(sqr(acc.m_v_sum_int)) / acc.m_v_n_int) << std::endl
		  << "External Variance: " << (acc.m_v_n_ext == 0 ? 0 : sum(acc.m_v_sum_ext_sqr) - sum(sqr(acc.m_v_sum_ext)) / acc.m_v_n_ext) << std::endl
		  << "Energy: " << acc.energy() << std::endl;
      }
    };

    /// \pre x and y must be canonical elements
    template <typename V>
    unsigned
    common_ancestor(unsigned x, unsigned y, const image2d< energy_t<V> >& e, const image2d<unsigned>& parent)
    {
      while (x != y) {
	if (e[x].depth > e[y].depth)
	  x = parent[x];
	else
	  y = parent[y];
      }
      return x;
    }

    template <typename Iterator, typename V, typename T>
    unsigned
    common_ancestor(Iterator x, const image2d< energy_t<V> >& e, const image2d<unsigned>& parent, const image2d<T>& K)
    {
      std::vector<unsigned> v;
      int mindepth = value_traits<int>::max();
      unsigned minp = 0;
      mln_forall(x) {
	if (!K.domain().has(x->point()))
	  continue;

	unsigned i = x->index();
	i = (parent[i] != (unsigned)-1 and K[i] == K[parent[i]]) ? parent[i] : i;
	v.push_back(i);
	if (e[i].depth < mindepth) {
	  mindepth = e[i].depth;
	  minp = i;
	}
      }

      bool modif = true;
      while (modif)
	{
	  modif = false;
	  for(unsigned& i: v)
	    if (e[i].depth > mindepth) {
	      i = parent[i];
	      modif = true;
	    } else if (e[i].depth == mindepth and i != minp) {
	      i = parent[i];
	      mindepth--;
	      minp = i;
	      modif = true;
	    }
	}
      mln_assertion(std::all_of(v.begin(), v.end(), [minp](unsigned x) { return x == minp; }));
      return v[0];
    }

  } // end of namespace internal


  template <typename V>
  auto
  norm(const V& v) -> typename std::enable_if< std::is_scalar<V>::value, decltype( std::abs(v) ) >::type
  {
    return std::abs(v);
  }

  template <typename V>
  auto
  norm(const V& v) -> typename std::enable_if< !std::is_scalar<V>::value, decltype(std::abs(v[0]) + std::abs(v[0])) >::type
  {
    auto s = std::abs(v[0]);
    for (int i = 1; i < V::ndim; ++i)
      s += std::abs(v[i]);
    return s;
  }



  template <typename V>
  image2d<float>
  curvature_on_edge(const image2d<V>& ima)
  {
    typedef decltype(V() + V()) Vec;
    typedef typename std::make_signed<Vec>::type vec_t;
    //static_assert(std::is_same<vec_t, int>::value, "here!");

    trace::entering("mln::curvature_on_edge");


    //auto norm = [] (const Vec& x) -> int { return std::abs(x[0] + x[1] + x[2]); };
    //auto norm = [] (const Vec& x) -> int { return std::abs(x); };
    auto sqr = [] (float x) -> float { return x*x; };

    auto domain = ima.domain();
    domain.pmin *= 2;
    domain.pmax = domain.pmax * 2 - 1;
    image2d<float> curv(domain, ima.border(), 0);



    mln_foreach(const point2d& p, ima.domain())
      {
	// Right edge
	{
	  float ux = norm( ((vec_t) ima.at(p + point2d{0,1}) - (vec_t) ima.at(p)));
	  float uy = norm( ((vec_t) ima.at(p + point2d{1,0}) - (vec_t) ima.at(p + point2d{-1,0}) +
			  (vec_t) ima.at(p + point2d{1,1}) - (vec_t) ima.at(p + point2d{-1,1})) / 4.0 );

	  //std::cout << p << " " << ux << " " << (- (vec_t) ima.at(p + point2d{-1,0})) << ":" << uy << std::endl;

	  float uxx = norm( (vec_t) ima.at(p + point2d{0,-1}) - (vec_t) ima.at(p)
			    -(vec_t) ima.at(p + point2d{0,1}) + (vec_t) ima.at(p + point2d{0,2})) / 2.0;

	  float uyy = norm( (vec_t) ima.at(p + point2d{-1,0}) + (vec_t) ima.at(p + point2d{-1,1})
			    -2 * ((vec_t) ima.at(p) + (vec_t) ima.at(p + point2d{0,1}))
			    + (vec_t) ima.at(p + point2d{1,0}) + (vec_t) ima.at(p + point2d{1,1}) ) / 2.0;

	  float uxy = norm( (vec_t) ima.at(p + point2d{1,0}) - (vec_t) ima.at(p + point2d{-1,0}) +
			    (vec_t) ima.at(p + point2d{1,1}) - (vec_t) ima.at(p + point2d{-1,1}) ) / 2.0;

	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{0,1};
	  if (den != 0)
	    curv.at(p_) = std::abs((uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den)));
	  else
	    curv.at(p_) = 0;

	  if (MLN_HAS_DEBUG and curv.at(p_) > 100 and p > point2d{1,1}) {
	    std::cout << p_ << std::endl;
	    io::imprint(ima | box2d{ p + point2d{-1,-1}, p + point2d{3,3}});

	    std::cout << norm((vec_t) ima.at(p + point2d{0,-1}) - (vec_t) ima.at(p)
			      -(vec_t) ima.at(p + point2d{0,1}) + (vec_t) ima.at(p + point2d{0,2})) << std::endl;

	    std::cout << curv.at(p_) << std::endl
		      << "ux: " << ux << std::endl
		      << "uy: " << uy << std::endl
		      << "uxx: " << uxx << std::endl
		      << "uyy: " << uyy << std::endl
		      << "uxy: " << uxy << std::endl
		      << "num: " << uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux) << std::endl
		      << "den: " << den * std::sqrt(den) << std::endl;
	    mln_assertion(false);
	  }
	}

	// Bottom edge
	{
	  float uy = norm( (vec_t) ima.at(p + point2d{1,0}) );

	  float ux = norm( (vec_t) ima.at(p + point2d{0,1}) - (vec_t) ima.at(p + point2d{0,-1}) +
			   (vec_t) ima.at(p + point2d{1,1}) - (vec_t) ima.at(p + point2d{1,-1}) ) / 4.0;

	  float uxx = norm( (vec_t) ima.at(p + point2d{-1,0}) - (vec_t) ima.at(p)
			    -(vec_t) ima.at(p + point2d{1,0}) + (vec_t) ima.at(p + point2d{2,0})) / 2.0;

	  float uyy = norm( (vec_t) ima.at(p + point2d{0,-1}) + (vec_t) ima.at(p + point2d{1,-1})
			    -2 * ((vec_t) ima.at(p) + (vec_t) ima.at(p + point2d{1,0}))
			    + (vec_t) ima.at(p + point2d{0,1}) + (vec_t) ima.at(p + point2d{1,1}) / 2.0);

	  float uxy = norm( (vec_t) ima.at(p + point2d{1,0}) - (vec_t) ima.at(p + point2d{0,-1}) +
			    (vec_t) ima.at(p + point2d{1,1}) - (vec_t) ima.at(p + point2d{1,-1})) / 2.0;

	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{1,0};
	  if (den != 0)
	    curv.at(p_) = std::abs((uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den)));
	  else
	    curv.at(p_) = 0;
	  // std::cout << curv.at(p_) << std::endl;
	  // assert(curv.at(p_) <= 1);
	}
      }

    trace::exiting();
    mln_postcondition(all(curv >= 0));
    //io::imsave(curv, "/tmp/curvature.tiff");
    return curv;
  }


  template <typename V, typename T>
  image2d<float>
  compute_energy_(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent_,
		  const std::vector<unsigned>& S, int eps = 5,
		  image2d< internal::energy_t<V> >* feedback = NULL)
  {
    typedef typename internal::energy_helper<V>::sum_type vec_sum_t;
    typedef typename internal::energy_helper<V>::sum_sqr_type vec_sum_sqr_t;

    extension::fill(ima, ima(ima.domain().pmin));

    image2d<float> curv = curvature_on_edge(ima);
    image2d<unsigned>& parent = const_cast<image2d<unsigned>&>(parent_);

    trace::entering("mln::compute_energy");

    image2d< internal::energy_t<V> > acc;
    resize(acc, K);

    //io::imprint(ima);
    //io::imprint(curv);

    // Compute depth attribute
    acc[S[0]].depth = 0;
    for (unsigned i = 1; i < S.size(); ++i)
      {
	unsigned x = S[i];
	if (K[x] != K[parent[x]]) // canonical element
	  acc[x].depth = acc[parent[x]].depth + 1;
      }


    // Compute attribute about edges
    {
      auto c4_idx = wrt_delta_index(curv, c4.dpoints);
      for (int i = S.size()-1; i > 0; --i)
	{
	  unsigned x = S[i];
	  if (K1::is_face_2(K.point_at_index(x))) {
	    acc[x].m_e_length += 4;
	    mln_foreach(int i, c4_idx) {
	      mln_assertion(K1::is_face_1(curv.point_at_index(x + i)));
	      acc[x].m_e_sumcurv += curv[x + i];
	    }
	  }
	  else if (K1::is_face_1(K.point_at_index(x))) {
	    acc[x].m_e_length += -2;
	    acc[x].m_e_sumcurv -= 2 * curv[x];
	  }
	  acc[parent[x]].m_e_length += acc[x].m_e_length;
	  acc[parent[x]].m_e_sumcurv += acc[x].m_e_sumcurv;
	  if (K[x] != K[parent[x]]) // real node assert !
	    assert(acc[x].m_e_sumcurv >= 0);
	}
      unsigned x = S[0];
      assert(K1::is_face_2(K.point_at_index(x)));
      acc[x].m_e_length += 4;
      mln_foreach(int i, c4_idx) {
	mln_assertion(K1::is_face_1(curv.point_at_index(x + i)));
	acc[x].m_e_sumcurv += curv[x + i];
      }
      // FIXME : Why the root has a negative curvature.
      acc[x].m_e_sumcurv = 0;
      //std::cout << acc[x].m_e_sumcurv << std::endl;
      //assert(acc[x].m_e_sumcurv >= 0);
    }



    typedef iterator_range< stditerator<std::vector<point2d>::const_iterator> > Vec;
    typedef dyn_neighborhood<std::vector<point2d>, dynamic_neighborhood_tag> Nbh;
    std::vector<point2d> dpoints;
    for (int i = -eps*2; i <= eps*2; i += 2)
      for (int j = -eps*2; j <= eps*2; j += 2)
	dpoints.emplace_back(i,j);

    Nbh ball(dpoints);

    // FIX TO HANDLE THE ROOT CORRECTLY
    parent[S[0]] = (unsigned) -1;

    // Compute attribute about internal area of components
    {
      mln_pixter(px, K);
      mln_iter(nx, ball(px));
      mln_forall(px)
      {
	if (K1::is_face_2(px->point()))
	  {
	    unsigned x = px->index();
	    x = (x != S[0] and K[x] == K[parent[x]]) ? parent[x] : x;

	    unsigned y = 0;
	    mln_forall(nx) {
	      if (!K.domain().has(nx->point()))
		{
		  y = (unsigned) -1; // This pixel is on the border
		  break;
		}
	    }

	    if (y == 0)
	      y = internal::common_ancestor(nx, acc, parent, K);

	    while (x != y)
	      {
		// p in Xint
		auto v = ima(px->point() / 2);
		acc[x].m_v_n_int += 1;
		acc[x].m_v_sum_int += (vec_sum_t) v;
		acc[x].m_v_sum_int_sqr += (vec_sum_sqr_t) v * (vec_sum_sqr_t) v;
		x = parent[x];
	      }
	  }
      }
    }

    // RESTORE ROOT
    parent[S[0]] = S[0];

    // Compute attribute about external area of components
    {
      mln_pixter(px, K);
      mln_iter(nx, ball(px));

      std::vector< std::pair<unsigned, unsigned> > branches;
      mln_forall(px)
      {
	if (K1::is_face_2(px->point()))
	  {
	    unsigned r = px->index();
	    r = (K[r] == K[parent[r]]) ? parent[r] : r;
	    branches.clear();

	    mln_forall(nx)
	    {
	      if (!K.domain().has(nx->point()))
		continue;

	      unsigned x = nx->index();
	      x = (K[x] == K[parent[x]]) ? parent[x] : x;

	      unsigned y = internal::common_ancestor(x, r, acc, parent);
	      branches.emplace_back(x,y);

	      while (x != y)
		{
		  if (acc[x].dejavu) {
		    x = parent[x];
		    continue;
		  }
		  // p in Xint
		  auto v = ima(px->point() / 2);
		  acc[x].m_v_n_ext += 1;
		  acc[x].m_v_sum_ext += (vec_sum_t) v;
		  acc[x].m_v_sum_ext_sqr += (vec_sum_sqr_t) v * (vec_sum_sqr_t) v;
		  acc[x].dejavu = true;
		  x = parent[x];
		}
	    }

	    unsigned x,y;
	    for (auto v: branches) {
	      std::tie(x,y) = v;
	      while (x != y) {
		acc[x].dejavu = false;
		x = parent[x];
	      }
	    }

	  }
      }
    }


    // Display
    // {
    //   image2d<bool> bdebug;
    //   resize(bdebug, K).init(false);

    //   for (int i = S.size()-1; i >= 0; --i)
    // 	{
    // 	  unsigned p = S[i];
    // 	  bdebug[p] = true;
    // 	  if (p == parent[p] or K[p] != K[parent[p]]) {
    // 	    std::cout << "=================" << std::endl
    // 		      << "Process: " << p << std::endl << acc[p] << std::endl;
    // 	    Kdisplay(bdebug);
    // 	  }
    // 	}
    // }


    // make enery image
    image2d<float> energy;
    resize(energy, acc);
    {
      energy[S[0]] = acc[S[0]].energy();
      for (unsigned x : S)
	if (K[x] != K[parent[x]])
	  energy[x] = acc[x].energy();
	else
	  energy[x] = energy[parent[x]];
    }

    trace::exiting();

    if (feedback != NULL)
      *feedback = acc;

    return energy;
  }

  template <typename V, typename T>
  image2d<float>
  compute_energy(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent_,
		  const std::vector<unsigned>& S, int eps)
  {
    return compute_energy_(ima, K, parent_, S, eps);
  }

  template <typename V, typename T>
  image2d<float>
  compute_energy(const image2d<V>& ima, const image2d<T>& K, const image2d<unsigned>& parent_,
		 const std::vector<unsigned>& S, image2d< internal::energy_t<V> > & acc, int eps = 5)
  {
    return compute_energy_(ima, K, parent_, S, eps, &acc);
  }


  template <typename V>
  image2d<float>
  close(const image2d<float>& energy, const image2d<V>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    image2d<float> imdilate = clone(energy);

    // Dilate
    for (int i = S.size()-1; i >= 0; --i)
      {
	unsigned x = S[i];
	unsigned q = parent[x];
	if (K[x] != K[q])
	  {
	    imdilate[x] = std::max(imdilate[x], energy[q]);
	    imdilate[q] = std::max(imdilate[q], energy[x]);
	  }
      }

    // Erode
    image2d<float> imerode = clone(imdilate);
    for (int i = S.size()-1; i >= 0; --i)
      {
	unsigned x = S[i];
	unsigned q = parent[x];
	if (K[x] != K[q])
	  {
	    imerode[x] = std::min(imerode[x], imdilate[q]);
	    imerode[q] = std::min(imerode[q], imdilate[x]);
	  }
      }

    return imerode;
  }

}

#endif // ! MUMFORD_SHAH_HPP
