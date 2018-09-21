#ifndef COMPUTE_ATTRIBUTE_HPP
# define COMPUTE_ATTRIBUTE_HPP

# include <mln/accu/accumulator.hpp>
# include <mln/core/algorithm/fill.hpp>
# include <mln/core/extension/fill.hpp>

namespace mln
{
  namespace accu
  {

    struct NoOpAccumulator
    {
      void init() {};

      template <typename T>
      void take(T) {}

      void take(const NoOpAccumulator&) {}

    };

  }

  namespace morpho
  {


    template <typename V, typename T, typename Accumulator2F = accu::NoOpAccumulator, typename Accumulator1F = accu::NoOpAccumulator>
    std::pair< image2d<Accumulator2F>, image2d<Accumulator1F> >
    tos_compute_attribute(const image2d<V>& K, const image2d<unsigned>& parent,
			  const std::vector<unsigned>& S, const image2d<T>& ori,
			  Accumulator2F acc2f = Accumulator2F (), Accumulator1F acc1f = Accumulator1F ())
    {
      static constexpr bool use_acc1f = !std::is_same<Accumulator1F, accu::NoOpAccumulator>::value;
      static constexpr bool use_acc2f = !std::is_same<Accumulator2F, accu::NoOpAccumulator>::value;

      static constexpr unsigned UNACTIVE = value_traits<unsigned>::max();

      mln_precondition(K.domain() == ori.domain());
      mln_precondition(K.domain() == parent.domain());

      struct edge_attribute
      {
	unsigned appear;
	unsigned vanish;
	float value;
      };

      box2d domain = K.domain();
      domain.pmax = domain.pmax * 2 - 1;
      image2d<edge_attribute> e;

      edge_attribute default_value = { UNACTIVE, UNACTIVE, 0.0};
      e.resize(domain, K.border(), default_value);
      extension::fill(e, default_value);

      // Kge of accumulators
      image2d<Accumulator2F> im_acc2f;
      image2d<Accumulator1F> im_acc1f;

      resize(im_acc2f, K).init(acc2f);
      resize(im_acc1f, K).init(acc1f);



      // Compute appear/vanish for edges
      // and 2-face attributes
      for (int i = S.size() - 1; i >= 0; --i)
	{
	  unsigned p  = S[i];
	  unsigned rp = K[parent[p]] == K[p] ? parent[p] : p;

	  point2d p_ = K.point_at_index(p);
	  point2d p2_ = 2 * p_;

	  if (K1::is_face_2(p_))
	    {
	      im_acc2f[rp].take(ori[p]);

	      mln_iter(n,  c4(p_));
	      mln_iter(n2, c4(p2_));
	      mln_forall(n, n2)
	      {
		if (e.at(*n2).appear == UNACTIVE) {
		  e.at(*n2).appear = rp;
		  vec3i grad = (vec3i) ori.at(*n) - (vec3i) ori[p];
		  e.at(*n2).value = (grad*grad).sum();
		} else
		  e.at(*n2).vanish = rp;
	      }
	    }
	  else if (K1::is_face_1v(p_))
	    {
	      mln_iter(n,  c2_h(p_));
	      mln_iter(n2, c2_h(p2_));
	      mln_forall(n, n2)
	      {
		//std::cout << "fuck" << p_ << *n << *n2 << std::endl;
		if (e.at(*n2).appear == UNACTIVE) {
		  e.at(*n2).appear = rp;
		  vec3i grad = (vec3i) ori.at(*n) - (vec3i) ori[p];
		  e.at(*n2).value = (grad*grad).sum();
		} else
		  e.at(*n2).vanish = rp;
	      }
	    }
	  else if (K1::is_face_1h(p_))
	    {
	      mln_iter(n,  c2_v(p_));
	      mln_iter(n2, c2_v(p2_));
	      mln_forall(n, n2)
	      {
		if (e.at(*n2).appear == UNACTIVE) {
		  e.at(*n2).appear = rp;
		  vec3i grad = (vec3i) ori.at(*n) - (vec3i) ori[p];
		  e.at(*n2).value = (grad*grad).sum();
		} else
		  e.at(*n2).vanish = rp;
	      }
	    }
	}

      //io::imprint_with_border(transform(e, [](edge_attribute x) { return x.appear; }));
      //io::imprint_with_border(transform(e, [](edge_attribute x) { return x.vanish; }));
      // Transmit edge atribute to nodes
      if (use_acc1f)
      {
	point2d pmin = e.domain().pmin;
	point2d pmax = e.domain().pmax;
	for (short i = pmin[0]-1; i < pmax[0]+1; ++i)
	  {
	    if (i % 4 == 2)
	      continue;
	    int step = i % 2 == 0 ? 2 : 4;
            for (short j = pmin[1] - (i % 2 == 0); j < pmax[1]+1; j += step) // foreach edges
	    {
	      point2d p{i,j};
	      unsigned cur_node  = e.at(p).appear;
	      unsigned last_node = e.at(p).vanish;
	      float    val = e.at(p).value;
	      assert(cur_node != UNACTIVE);
	      if (last_node == UNACTIVE)
		im_acc1f[cur_node].take(val);
	      else
		while (cur_node != last_node) {
		  im_acc1f[cur_node].take(val);
		  cur_node = parent[cur_node];
		}
	    }
	  }
      }

      if (true)
	{
	  for (unsigned p: S)
	    {
	      unsigned q = parent[p];
	      if (K[p] == K[q])
		{
		  if (use_acc2f)
		    im_acc2f[p] = im_acc2f[q];
		  if (use_acc1f)
		    im_acc1f[p] = im_acc1f[q];
		}
	    }
	}

      return std::make_pair(im_acc2f, im_acc1f);
    }

  }

}

# ifndef MLN_INCLUDE_ONLY

# endif // ! MLN_INCLUDE_ONLY

#endif // ! COMPUTE_ATTRIBUTE_HPP
