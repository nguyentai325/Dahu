#ifndef MLN_CORE_MORPHO_MAXTREE_PQUEUE_HPP
# define MLN_CORE_MORPHO_MAXTREE_PQUEUE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/morpho/pqueue_fast.hpp>
# include <vector>
# include <stack>
# include <queue>


#include <mln/io/imprint.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {

      template <typename V, typename Neighborhood, typename StrictWeakOrdering, bool parallel>
      void
      maxtree_flood_pqueue_algorithm(const image2d<V>& ima,
				     image2d<typename image2d<V>::size_type>& parent,
				     const Neighborhood& , StrictWeakOrdering cmp, typename image2d<V>::size_type* Send)
      {
	typedef typename image2d<V>::size_type size_type;
	static constexpr bool use_dejavu = true;
	static constexpr size_type UNINITIALIZED = std::numeric_limits<size_type>::max();
	static constexpr size_type INQUEUE = 0;
	static constexpr std::size_t nlevels = (std::size_t)1 << value_traits<V>::quant;


	image2d<bool> deja_vu;
	typedef std::vector<size_type> vec_t;
	vec_t v_stack; v_stack.reserve(std::min<std::size_t>(ima.domain().size(), nlevels));
	//vec_t v_pqueue; v_pqueue.reserve(ima.domain().size());
	//auto fcmp = [&ima, cmp] (size_type x, size_type y) { return cmp(ima[x], ima[y]); };
	std::stack<size_type, vec_t> stack(std::move(v_stack));
	priority_queue_ima<V, StrictWeakOrdering> pqueue(ima, cmp);

	  //std::priority_queue<size_type, vec_t, decltype(fcmp)>
	  //pqueue(fcmp, std::move(v_pqueue));
	auto nbh_delta_indexes = wrt_delta_index(ima, Neighborhood::dpoints);

	size_type* Send_ = Send;
	(void) Send_;

	size_type first_index = ima.index_of_point(ima.domain().pmin);
	size_type last_index = ima.index_of_point(ima.domain().pmax) - ima.index_strides()[0];


	// INIT
	{
	  if (use_dejavu) {
	    resize(deja_vu, ima).init(false);
	    extension::fill(deja_vu, true);
	  }

	  // Get min element and reserve queue
	  size_type pmin = ima.index_of_point(ima.domain().pmin);
	  //V vmin = ima[pmin];
	  {
	    // mln_pixter(px, ima);
	    // mln_forall(px)
	    // {
	    //   if (cmp(px->val(), vmin))
	    // 	{
	    // 	  vmin = px->val();
	    // 	  pmin = px->index();
	    // 	}
	    // }

	    pqueue.push(pmin);
	    stack.push(pmin);
	    if (use_dejavu)
	      deja_vu[pmin] = true;
	    else
	      parent[pmin] = 0;
	  }
	}


	while (!pqueue.empty())
	  {
	  flood:
	    size_type p = pqueue.top();
	    size_type repr = stack.top();
	    assert(ima[p] == ima[repr]);

	    mln_foreach(auto k, nbh_delta_indexes)
	      {
		auto q = p + k;
		bool processed;
		if (use_dejavu)
		  processed = deja_vu[q];
		else if (!parallel)
		  processed = (parent[q] != UNINITIALIZED);
		else
		  processed = !(first_index <= q and q < last_index) or (parent[q] != UNINITIALIZED); 

		if (!processed) {
		  pqueue.push(q);
		  if (use_dejavu)
		    deja_vu[q] = true;
		  else
		    parent[q] = INQUEUE;
		  if (cmp(ima[p], ima[q])) {
		    stack.push(q);
		    goto flood;
		  }
		}
	      }

	    // p done
	    pqueue.pop();
	    parent[p] = repr;
	    if (!parallel and p != repr) *(--Send) = p;

	    // Stack handle : attach to parent
	    if (pqueue.empty()) break;
	    size_type next = pqueue.top();
	    assert(!cmp(ima[repr], ima[next]));
	    if (cmp(ima[next], ima[repr]))
	      {
		//size_type par;
		stack.pop();
		while (!stack.empty() and cmp(ima[next], ima[stack.top()]))
		  {
		    if (!parallel) *(--Send) = repr;
		    repr = (parent[repr] = stack.top());
		    stack.pop();
		  }
		if (stack.empty() or cmp(ima[stack.top()], ima[next]))
		  stack.push(next);

		parent[repr] = stack.top();
		if (!parallel) *(--Send) = repr;
	      }
	  }

	assert(stack.size() == 1);
	if (!parallel)
	  *(--Send) = stack.top();

	// std::cout << ima.domain().size()  << std::endl;
	// std::cout << cpt  << std::endl;
	// std::cout << cpt2  << std::endl;
	// std::cout << "D:" << (cpt2-cpt)  << std::endl;
	// std::cout << cpt3  << std::endl;
	// std::cout << cpt4  << std::endl;

	if (!parallel)
	  assert((Send + ima.domain().size()) == Send_);


	/* 	// Flood
	while (!pqueue.empty())
	  {
	    size_type p = pqueue.top();
	    size_type repr = stack.top();

	    if (cmp(ima[p], ima[repr]))
	      {
		// the component represented by repr is over
		// we attach repr to its parent
		stack.pop();
		assert(!stack.empty());
		size_type par = stack.top();
		while (cmp(ima[p], ima[par]))
		  {
		    stack.pop();
		    parent[repr] = par;
		    if (!parallel) *(--Send) = repr;
		    repr = par;
		    par = stack.top();
		  }
		if (cmp(ima[par], ima[p])) { // ima[par] < ima[p]
		  stack.push(p);
		  par = p;
		}
		parent[repr] = par;
		if (!parallel) *(--Send) = repr;
		repr = par;
	      }


	    assert(ima[repr] == ima[p]);
	    bool done = true;
	    mln_foreach(auto k, nbh_delta_indexes)
	      {
		auto q = p + k;
		bool processed;
		if (use_dejavu)
		  processed = deja_vu[q];
		else if (!parallel)
		  processed = (parent[q] != UNINITIALIZED);
		else
		  processed = !(first_index <= q and q < last_index) or (parent[q] != UNINITIALIZED); 

		if (!processed) {
		  pqueue.push(q);
		  if (use_dejavu)
		    deja_vu[q] = true;
		  else
		    parent[q] = INQUEUE;
		  if (cmp(ima[p], ima[q])) {
		    stack.push(q);
		    done = false;
		    break;
		  }
		}
	      }

	    if (done) {
	      parent[p] = repr;
	      if (!parallel and p != repr) *(--Send) = p;
	      pqueue.pop();
	    }
	  }

	if (!parallel)
	  while (!stack.empty()) {
	    *(--Send) = stack.top();
	    stack.pop();
	  }

	if (!parallel)
	  assert((Send + ima.domain().size()) == Send_);
	*/
      }

    } // end of namespace mln::morpho::internal


    template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
    std::pair<image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
    maxtree_pqueue(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
    {
      typedef typename image2d<V>::size_type size_type;
      image2d<size_type> parent;
      resize(parent, ima);

      std::vector<size_type> S(ima.domain().size());

      internal::maxtree_flood_pqueue_algorithm(ima, parent, nbh, cmp, &(*S.rend));
      return parent;
    }

  }  // end of namespace mln::morpho

} // end of namespace mln

#endif // ! MLN_CORE_MORPHO_MAXTREE_PQUEUE_HPP
