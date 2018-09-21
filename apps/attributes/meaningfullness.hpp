#ifndef APPS_ATTRIBUTES_MEANINGFULLNESS_HPP
# define MUMFORD_SHAH_HPP

# include <mln/core/trace.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/core/neighb2d.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/wrt_offset.hpp>

# include <mln/core/colors.hpp>
# include <mln/core/math_ops.hpp>
# include <mln/core/vec/vec_math_ops.hpp>
# include <apps/tos/topology.hpp>

namespace mln
{


  /// \brief Compute the meaningfullness energy for each nodes.
  /// \pre K must be twice as big as ima.
  ///
  /// This method computes the energy defined in xu.2012.icip
  ///
  /// E = α.Eᵢ + Eₑ + β.Eₓ
  /// where:
  /// * Eᵢ is the internal energy defined as
  ///   Eᵢ(A) = (V_i(A,ε) + V\_{ext}(A,ε)) / V(A,ε)
  ///
  ///
  template <typename V, typename T>
  image2d<float>
  meaningfullness(const image2d<V>& ima,
		  const image2d<T>& K,
		  const image2d<unsigned>& parent,
		  const std::vector<unsigned>& S,
		  float alpha,
		  float beta,
		  float gamma,
		  int eps = 5);


  /// \brief Take an image and puts its curvature on the edges
  /// The image is doubled, non-edges pixels have undefined values.
  /// If the input value type is a rgb, the image is first converted
  /// to gray-levels (by averaging r,g,b values)
  template <typename V>
  image2d<float>
  curvature_on_edge(const image2d<V>& ima);


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

    template <typename V, bool scalar = std::is_scalar<V>::value >
    struct energy_helper
    {
      enum { dim = 1 };
      typedef uint64 sum_type;
      typedef uint64 sum_sqr_type;
    };

    template <typename V>
    struct energy_helper<V, false>
    {
      enum { dim = V::ndim };
      typedef vec<uint64, V::ndim> sum_type;
      typedef vec<uint64, V::ndim> sum_sqr_type;
    };

    template <typename V>
    struct energy_t
    {
      typedef typename energy_helper<V>::sum_type sum_type;
      typedef typename energy_helper<V>::sum_sqr_type sum_sqr_type;
      enum { dim = energy_helper<V>::dim };

      static float alpha;  //  Static variable ? really ?
      static float beta;
      static float gamma;

      bool		dejavu = false;
      int		depth = -1;
      int		m_e_length = 0;
      float		m_e_sumcurv = 0;
      int		m_v_n_int = 0;
      sum_type		m_v_sum_int = literal::zero;
      sum_sqr_type	m_v_sum_int_sqr = literal::zero;
      int		m_v_n_ext = 0;
      sum_type		m_v_sum_ext = literal::zero;
      sum_sqr_type	m_v_sum_ext_sqr = literal::zero;
      unsigned		m_area = 0;

      energy_t()
      {
      }

      template <typename T>
      static float cvar(T x, T x2, float n)
      {
	if (n == 0)
	  return 0;
	return sum(x2 / n - sqr(x / n)) / (float) dim;
      }

      /// Return the variance σ² of the gaussian distribution
      /// with a covariance matrix: V = σ².Idⁿ where n is the dimension.
      static float var(sum_type x, sum_sqr_type x2, float n)
      {
	if (n == 0) return 0;

	// float s = sum(x2 / n) - sum(sqr(x / n));
	// which is equivalent to
	float s = sum(x2 / n - sqr(x) / (n*n)  );
	return s / (float) dim; //std::pow(s, 1.0 / dim);
      }


      static
      std::pair<float, float>
      inter(float m1, float var1,
	    float m2, float var2, float& delta)
      {
	float a = var2 - var1;
	float b = -2 * (m1 * var2 - m2 * var1);
	float c = sqr(m1)*var2 - sqr(m2)*var1 - 2 * var1 * var2 * std::log(sqrt(var2/var1));
	delta = sqr(b) - 4*a*c;
	if (delta == 0)
	  return std::make_pair(-b / (2*a), 0.0f);
	float r1 = (-b - sqrt(delta)) / (2*a);
	float r2 = (-b + sqrt(delta)) / (2*a);
	return std::make_pair(r1, r2);
      }


      float external_energy() const {
	float Vin = var(m_v_sum_int, m_v_sum_int_sqr, m_v_n_int);
	float Vout = var(m_v_sum_ext, m_v_sum_ext_sqr, m_v_n_ext);
	// float Vtotal = var(m_v_sum_int + m_v_sum_ext,
	// 		   m_v_sum_int_sqr + m_v_sum_ext_sqr,
	// 		   m_v_n_int + m_v_n_ext);

	auto c1_ = (m_v_sum_int / (float)m_v_n_int);
	auto c2_ = (m_v_sum_ext / (float)m_v_n_ext);
	auto v_ = c2_ - c1_;
	auto v = (v_) / (float)l2norm(v_);

	float c1 = sum(c1_ * v);
	float c2 = sum(c2_ * v);

	//float dinter = sum(sqr(c1-c2));

	if (Vin == 0 or Vout == 0)
	  return 1.0;

	// KULLBACK-LEIBER HERE
	/*
	auto KL = [dinter] (float var1, float var2) {
	  float std1 = std::sqrt(var1);
	  float std2 = std::sqrt(var2);
	  return std::log(std2/std1) + (var1 + dinter) / (2 * var2) - 0.5;
	};
	float dKL = (KL(Vin, Vout) + KL(Vout, Vin)) * 0.5;
	*/

	// HELLINGER HERE
	/*
	float std1 = std::sqrt(Vin);
	float std2 = std::sqrt(Vout);
	float Hellinger = 1.0 - std::sqrt(2*(std1*std2)/(Vin + Vout)) * std::exp(-0.25 * dinter / (Vin + Vout));
	return 1.0 - Hellinger;
	*/

	// min / max integral proba HERE
	{
	  float delta;
	  float a,b;
	  auto hf = [c1, Vin] (float x) { return 0.5 * (1 + std::erf( (x-c1)/std::sqrt(2*M_PI*Vin))); };
	  auto hg = [c2, Vout] (float x) { return 0.5 * (1 + std::erf( (x-c2)/std::sqrt(2*M_PI*Vout))); };

	  std::tie(a,b) = inter(c1, Vin, c2, Vout, delta);
	  if (delta == 0)
	    {
	      float d = std::abs(hf(a) - hg(a));
	      return (1-d)/(1+d);
	    }
	  else
	    {
	      float if1, if2, if3, ig1, ig2, ig3;
	      if1 = hf(a);
	      ig1 = hg(a);
	      if2 = hf(b) - if1;
	      ig2 = hg(b) - ig1;
	      if3 = 1 - if1 - if2;
	      ig3 = 1 - ig1 - ig2;
	      float m1, M1, m2, M2, m3, M3;
	      std::tie(m1, M1) = std::minmax(if1, ig1);
	      std::tie(m2, M2) = std::minmax(if2, ig2);
	      std::tie(m3, M3) = std::minmax(if3, ig3);
	      return (m1+m2+m3)/(M1+M2+M3);
	    }
	}


	// return min(Vin, Vout) / Vtotal;
	// float Vtotal = cvar(m_v_sum_int + m_v_sum_ext,
	// 		    m_v_sum_int_sqr + m_v_sum_ext_sqr,
	// 		    m_v_n_int + m_v_n_ext);
	// if (Vtotal == 0)
	//   return 0;
	//return (Vin + Vout) / Vtotal;
      }

      float energy() const
      {
	//auto sqr = [](float x) { return x*x; };
	float Eext = external_energy();
	float Eint = m_e_sumcurv / m_e_length;
	float Econ = std::exp(-gamma * m_area);
	// std::cout << "Contour length: " <<  m_e_length << std::endl
	// 	  << "Curvature: " << m_e_sumcurv << std::endl;
	assert(m_e_length > 0);
	assert(Eext >= 0);
	assert(Eint >= 0);
	assert(Econ >= 0);
	return std::max(alpha * Eint + Eext, beta * Econ);
      }

      float energy(float alpha_, float beta_, float gamma_ = 1.0) const
      {
	//auto sqr = [](float x) { return x*x; };
	float Eext = external_energy();
	float Eint = m_e_sumcurv / m_e_length;
	float Econ = std::exp(- gamma_ * m_area);
	// std::cout << "Contour length: " <<  m_e_length << std::endl
	// 	  << "Curvature: " << m_e_sumcurv << std::endl;
	assert(m_e_length > 0);
	assert(Eext >= 0);
	assert(Eint >= 0);
	assert(Econ >= 0);
	return std::max(alpha_ * Eint + Eext, beta_ * Econ);
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
		  << "Internal Variance: " << acc.cvar(acc.m_v_sum_int, acc.m_v_sum_int_sqr, acc.m_v_n_int) << std::endl
		  << "External Variance: " << acc.cvar(acc.m_v_sum_ext, acc.m_v_sum_ext_sqr, acc.m_v_n_ext) << std::endl
		  << "Energy: " << acc.energy() << std::endl;
      }
    };

    template <typename T>
    float energy_t<T>::alpha;

    template <typename T>
    float energy_t<T>::beta;

    template <typename T>
    float energy_t<T>::gamma;

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
  curvature_on_edge(const image2d<V>& ima_)
  {
    static_assert(std::is_same<V, uint8>::value or
		  std::is_same<V, rgb8>::value,
		  "V must either uint8 or rgb8.");

    struct to_uint8 {
      static
      image2d<uint8>
      foo(const image2d<uint8>& x)
      {
	return x;
      }

      static
      image2d<uint8>
      foo(const image2d<rgb8>& ima_)
      {
	return transform(ima_, [] (rgb8 x) -> uint8 { return (x[0] + x[1] + x[2]) / 3; });
      }
    };

    image2d<uint8> ima = to_uint8::foo(ima_);


    trace::entering("mln::curvature_on_edge");

    auto domain = ima.domain();
    domain.pmin *= 2;
    domain.pmax = domain.pmax * 2 - 1;
    image2d<float> curv(domain, ima.border(), 0);


    //auto l1norm = [](vec_t x) { return x; };

    mln_foreach(const point2d& p, ima.domain())
      {
	// Right edge
	{
	  float ux = ima.at(p + point2d{0,1}) - ima.at(p);
	  float uy =
	    (ima.at(p + point2d{1,0}) - ima.at(p + point2d{-1,0}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{-1,1})) / 4.0;

	  float uxx =
	    (ima.at(p + point2d{0,-1}) - ima.at(p) -
	     ima.at(p + point2d{0,1})  + ima.at(p + point2d{0,2})) / 2.0;

	  float uxy =
	    (-ima.at(p + point2d{1,0}) + ima.at(p + point2d{-1,0}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{-1,1})) / 2.0;


	  float uyy =
	    (ima.at(p + point2d{-1,0}) + ima.at(p + point2d{-1,1}) +
	     ima.at(p + point2d{1,0}) + ima.at(p + point2d{1,1})
	     - 2 * ima.at(p) - 2 * ima.at(p + point2d{0,1})) / 2.0;


	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{0,1};
	  if (den != 0)
	    curv.at(p_) = std::abs(uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den));
	  else
	    curv.at(p_) = 0;


	  // if (MLN_HAS_DEBUG and curv.at(p_) > 100 and p > point2d{1,1}) {
	  //   std::cout << p_ << std::endl;
	  //   io::imprint(ima | box2d{ p + point2d{-1,-1}, p + point2d{3,3}});


	  //   std::cout << curv.at(p_) << std::endl
	  // 	      << "ux: " << ux << std::endl
	  // 	      << "uy: " << uy << std::endl
	  // 	      << "uxx: " << uxx << std::endl
	  // 	      << "uyy: " << uyy << std::endl
	  // 	      << "uxy: " << uxy << std::endl
	  // 	      << "num: " << uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux) << std::endl
	  // 	      << "den: " << den * std::sqrt(den) << std::endl;
	  //   mln_assertion(false);
	  // }
	}

	// Bottom edge
	{
	  float uy = ima.at(p + point2d{1,0}) - ima.at(p);

	  float ux =
	    (ima.at(p + point2d{0,1}) - ima.at(p + point2d{0,-1}) +
	     ima.at(p + point2d{1,1}) - ima.at(p + point2d{1,-1})) / 4.0;

	  float uyy =
	    (ima.at(p + point2d{-1,0}) - ima.at(p) -
	     ima.at(p + point2d{1,0}) + ima.at(p + point2d{2,0})) / 2.0;

	  float uxx =
	    (ima.at(p + point2d{0,-1}) + ima.at(p + point2d{1,-1})
	     - 2 * ima.at(p) - 2 * ima.at(p + point2d{1,0}) +
	     ima.at(p + point2d{0,1}) + ima.at(p + point2d{1,1})) / 2.0;

	  float uxy =
	    (ima.at(p + point2d{0,-1}) - ima.at(p + point2d{0,1}) -
	     ima.at(p + point2d{1,-1}) + ima.at(p + point2d{1,1})) / 2.0;

	  float den = (sqr(ux) + sqr(uy));
	  point2d p_ = p * 2 + point2d{1,0};
	  if (den != 0)
	    curv.at(p_) = std::abs(uxx * sqr(uy) - 2 * uxy * ux *uy + uyy * sqr(ux)) / (den * std::sqrt(den));
	  else
	    curv.at(p_) = 0;
	}
      }

    trace::exiting();
    mln_postcondition(all(curv >= 0));


    //io::imsave(tmp, "/tmp/curvature.tiff");
    return curv;
  }


  template <typename V, typename T>
  image2d<float>
  compute_energy_(const image2d<V>& ima,
		  const image2d<T>& K,
		  const image2d<unsigned>& parent_,
		  const std::vector<unsigned>& S,
		  float alpha,
		  float beta,
		  float gamma,
		  int eps = 5,
		  image2d< internal::energy_t<V> >* feedback = NULL)
  {
    typedef typename internal::energy_helper<V>::sum_type vec_sum_t;
    typedef typename internal::energy_helper<V>::sum_sqr_type vec_sum_sqr_t;

    extension::fill(ima, ima(ima.domain().pmin));

    image2d<float> curv = curvature_on_edge(ima);
    image2d<unsigned>& parent = const_cast<image2d<unsigned>&>(parent_);

    trace::entering("mln::compute_energy");

    internal::energy_t<V>::alpha = alpha;
    internal::energy_t<V>::beta = beta;
    internal::energy_t<V>::gamma = gamma;

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


    mln_assertion( are_indexes_compatible(curv, parent) );
    mln_assertion( are_indexes_compatible(curv, K) );


    // Compute attribute about edges
    {
      auto c4_idx = wrt_delta_index(curv, c4.dpoints);
      for (int i = S.size()-1; i > 0; --i)
	{
	  unsigned x = S[i];
	  if (K1::is_face_2(K.point_at_index(x))) {
	    acc[x].m_area += 1;
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
	  acc[parent[x]].m_area += acc[x].m_area;
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
      energy[S[0]] = acc[S[0]].energy(alpha, beta, gamma);
      for (unsigned x : S)
	if (K[x] != K[parent[x]])
	  energy[x] = acc[x].energy(alpha, beta, gamma);
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
  meaningfullness(const image2d<V>& ima,
		  const image2d<T>& K,
		  const image2d<unsigned>& parent_,
		  const std::vector<unsigned>& S,
		  float alpha,
		  float beta,
		  float gamma,
		  int eps)
  {
    return compute_energy_(ima, K, parent_, S, alpha, beta, gamma, eps);
  }

  template <typename V, typename T>
  image2d<float>
  meaningfullness(const image2d<V>& ima,
		  const image2d<T>& K,
		  const image2d<unsigned>& parent_,
		  const std::vector<unsigned>& S,
		  image2d< internal::energy_t<V> >& acc,
		  float alpha,
		  float beta,
		  float gamma,
		  int eps = 5)
  {
    return compute_energy_(ima, K, parent_, S, alpha, beta, gamma, eps, &acc);
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
