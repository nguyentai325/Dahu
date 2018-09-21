#ifndef MLN_MORPHO_MAXTREE_MAXTREE_QUEUE_HPP
# define MLN_MORPHO_MAXTREE_MAXTREE_QUEUE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/fill.hpp>

# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/pqueue_fast.hpp>

# include <mln/core/wrt_offset.hpp>

# include <vector>
# include <stack>
# include <queue>

namespace mln
{

  namespace morpho
  {

    namespace impl
    {

      template <typename I, typename Im, typename Pa, typename Neighborhood, typename StrictWeakOrdering>
      component_tree<typename I::size_type, mln_ch_value(I, unsigned)>
      maxtree_queue_indexes(const I& ima, const Im& ima1, const Pa& ima2, const Neighborhood& nbh, StrictWeakOrdering cmp)
      {
	typedef typename I::size_type size_type;
	typedef mln_value(I) V;


	// 1. Create the component tree
	// Allocate enough memory to prevent reallocation
	// {
	typedef mln_ch_value(I, unsigned) map_t;
	typedef morpho::internal::component_tree_node node_t;
	typedef component_tree<size_type, map_t> tree_t;

	component_tree<size_type, mln_ch_value(I, size_type)> ctree;

	auto& nodes = ctree._get_data()->m_nodes;
	auto& S     = ctree._get_data()->m_S;
	auto& pmap  = ctree._get_data()->m_pmap;
	auto& Uv  = ctree._get_data()->m_Uv;
	auto& parent_pixel  = ctree._get_data()->m_parent_pixel;

	size_type sz = ima.domain().size();
	Uv = ima1;
	parent_pixel = ima2;


	std::cout <<  "sz   "   << sz   << std::endl;
	
	nodes.resize(sz + 1);   // grow from the back
	S.resize(sz);		    //
	pmap = imchvalue<size_type>(ima);	//


	size_type spos = sz;
	// }

	// 1.b some methods to handle the nodes array like a stack
	// {
	size_type stack_top_position = 0;
	size_type npos = 1; // The current nodes vector size

	auto stack_empty = [&stack_top_position] () {
	  return stack_top_position == 0;
	};
	auto stack_top = [&nodes, &stack_top_position] () -> node_t& {
	  mln_precondition(stack_top_position != 0);
	  return nodes[stack_top_position];
	};
	auto stack_push = [&nodes, &stack_top_position, &npos] (const node_t& x) {
	  mln_assertion(stack_top_position < npos);
	  nodes[npos] = x;
	  nodes[npos].m_parent = stack_top_position;
	  stack_top_position = npos++;
	};
	auto stack_pop = [&nodes, &stack_top_position] () {
	  mln_precondition(stack_top_position != 0);
	  stack_top_position = nodes[stack_top_position].m_parent;
	};

	// }


	// 2. Create auxiliary structure for the computation
	// {
	mln_ch_value(I, bool) deja_vu = imchvalue<bool>(ima).adjust(nbh).init(false);
	priority_queue_ima<V, StrictWeakOrdering> pqueue(ima, cmp);

	extension::fill(deja_vu, true);
	// }

	// 3. Initialization
	{
	  size_type pmin = ima.index_of_point(ima.domain().pmin);
	  pqueue.push(pmin);
	  stack_push( node_t{0,0,0, npos, pmin} );
	  deja_vu[pmin] = true;
	}


	// 4. algo
	// next_sibling is actually the last node in the substree
	// it will be corrected afterward

	// Fixme non-generic -> only dynamic neighborhood
	auto nbh_delta_indexes = wrt_delta_index(ima, Neighborhood::dpoints);


	while (!pqueue.empty())
	  {
	  flood:
	    size_type p = pqueue.top();
	    // size_type repr = stack.top();
	    // mln_assertion(not cmp(ima[p], ima[repr]) and not cmp(ima[repr], ima[p]));

	    V clevel = ima[p];

	    mln_foreach(auto k, nbh_delta_indexes)
	      {
		auto q = p + k;
		bool processed = deja_vu[q];

		if (!processed) {
		  pqueue.push(q);
		  deja_vu[q] = true;
		  if (cmp(clevel, ima[q])) {
		    stack_push( node_t{0,0,0,npos,q} );
		    goto flood;
		  }
		}
	      }

	    // p done
	    // insert p in S and set him as representative
	    // for the current node if it is note yet defined

	    pqueue.pop();
	    S[--spos] = p;
	    pmap[p] = stack_top_position;

	    // If this is the last point at this
	    // level: insert the node
	    if (pqueue.empty()) break; // We finished to process the queue

	    if (cmp(ima[pqueue.top()], clevel))
	      {
		size_type next = pqueue.top();
		size_type current = stack_top_position;

		stack_pop();
		nodes[current].m_point_index = spos;

		while (!stack_empty() and cmp(ima[next], ima[stack_top().m_point_index]))
		  {
		    size_type q = stack_top_position;
		    nodes[current].m_parent = q; //no-op
		    nodes[current].m_prev = nodes[q].m_next_sibling;
		    nodes[nodes[q].m_next_sibling].m_next = current;
		    nodes[q].m_next_sibling = nodes[current].m_next_sibling;

		    current = stack_top_position;
		    stack_pop();
		  }
		if (stack_empty() or cmp(ima[stack_top().m_point_index], ima[next]))
		  stack_push(node_t{0,0,0,npos, next});

		size_type q = stack_top_position;
		nodes[current].m_parent = q; //no-op
		nodes[current].m_prev = nodes[q].m_next_sibling;
		nodes[nodes[q].m_next_sibling].m_next = current;
		nodes[q].m_next_sibling = nodes[current].m_next_sibling;
	      }
	  }
	mln_assertion( stack_top().m_parent == 0);
	// Handle the root
	nodes[stack_top_position].m_point_index = 0;

	// Sentinel
	nodes[0] = node_t{0,  // parent -> itself
			  nodes[stack_top_position].m_next_sibling, // prev = last node
			  0,  // next -> itself
			  0,  // next_sibling -> itself
			  sz };

	for (unsigned i = 1; i < npos; ++i)
	  nodes[i].m_next_sibling = nodes[nodes[i].m_next_sibling].m_next;

	nodes.resize(npos);
	nodes.shrink_to_fit();

	return ctree.get_subtree(stack_top_position);
      }

    } // end of namespace mln::morpho::impl

  }

}

#endif // ! MLN_MORPHO_MAXTREE_MAXTREE_QUEUE_HPP
