#ifndef MLN_CORE_MORPHO_MAXTREE_HQUEUE_HPP
# define MLN_CORE_MORPHO_MAXTREE_HQUEUE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/value/indexer.hpp>
# include <mln/morpho/bounded_hqueue.hpp>


#include <mln/io/imprint.hpp>

namespace mln
{

  namespace morpho
  {

    namespace internal
    {

      template <typename I, typename Neighborhood, typename StrictWeakOrdering, bool parallel>
      struct maxtree_flood_algorithm
      {
	typedef mln_value(I) V;
	typedef typename I::size_type size_type;

	typedef typename indexer<V,StrictWeakOrdering>::index_type level_t;
	typedef size_type  elt_type; // either point or index
	static constexpr std::size_t nlevels = 1ul << value_traits<level_t>::quant;
	static constexpr bool use_dejavu = true;

	static constexpr size_type UNINITIALIZED =  std::numeric_limits<size_type>::max();
	static constexpr size_type INQUEUE = 0;

	level_t flood(level_t level)
	{
	  mln_precondition(m_has_repr[level]);
	  mln_precondition(not m_q.empty(level));
	  elt_type r = m_repr[level];
	  mln_precondition(m_h(m_ima[r]) == level);

          //std::cout << "Start of flood: " << (int)level << "(repr:" << r << ")" << std::endl;
          //io::imprint(m_parent);

	  while (!m_q.empty(level))
	    {
	      size_type p = m_q.pop_at_level(level);
	      if (!parallel and p != r) { *(--m_out) = p; }
	      m_parent[p] = r;

	      mln_foreach(auto k, m_nbh_delta_indexes)
		{
		  auto q = p + k;

		  bool processed;
		  if (use_dejavu)
		    processed = m_deja_vu[q];
		  else if (!parallel) // no bound checking
		    processed = (m_parent[q] != UNINITIALIZED);
		  else
		    processed = q < m_first_index or m_last_index <= q or (m_parent[q] != UNINITIALIZED);

		  if (!processed)
		    {
		      level_t newlevel = m_h(m_ima[q]);
		      if (!m_has_repr[newlevel])
			{
			  m_repr[newlevel] = q;
			  m_has_repr[newlevel] = true;
			}
                      //std::cout << "++ push: " << q << " @ " << (int)newlevel << std::endl;
		      m_q.push_at_level(q, newlevel);
		      if (use_dejavu)
			m_deja_vu[q] = true;
		      else
			m_parent[q] = INQUEUE;

		      if (level < newlevel)
			do {
			  newlevel = flood(newlevel);
			} while (level < newlevel);
		    }
		}
	    }

	  // Attach to parent
	  if (!parallel) { *(--m_out) = r; }
	  m_has_repr[level] = false;
	  while (level > value_traits<level_t>::min()) {
            --level;
            if (m_has_repr[level]) {
              m_parent[r] = m_repr[level];
              break;
            }
	  }

	  return level;
	}

	maxtree_flood_algorithm(const I& ima, image2d<size_type>& parent,
                                const Neighborhood&, StrictWeakOrdering, size_type* Send = NULL)
	  : m_ima(ima), m_parent (parent),
	    m_nbh_delta_indexes(wrt_delta_index(ima, Neighborhood::dpoints)),
	    m_has_repr {false,}, m_out (Send)
	{
	  if (use_dejavu) {
	    resize(m_deja_vu, ima).init(false);
	    extension::fill(m_deja_vu, true);
	  }

	  m_first_index = m_ima.index_of_point(m_ima.domain().pmin);
	  m_last_index =  m_ima.index_of_point(m_ima.domain().pmax) - ima.index_strides()[0];

	  // Get min element and reserve queue
	  size_t pmin = ima.index_of_point(ima.domain().pmin);
	  level_t vmin = value_traits<level_t>::max();
	  {
	    std::vector<std::size_t> h(nlevels, 0);
	    //size_type h[nlevels] = {0,};

	    mln_pixter(px, ima);
	    mln_forall(px)
	    {
	      level_t l = m_h(px->val());
	      ++h[l];
	      if (l < vmin)
		{
		  vmin = l;
		  pmin = px->index();
		}
	    }

	    m_q.init(&h[0]);
	  }

	  // Start flooding
	  //std::cout << (int) vmin << "@" << pmin << std::endl;
	  m_q.push_at_level(pmin, vmin);
	  m_repr[vmin] = pmin;
	  m_has_repr[vmin] = true;
	  if (use_dejavu)
	    m_deja_vu[pmin] = true;
	  else
	    m_parent[pmin] = INQUEUE;
	  flood(vmin);
	}

	static void run(const I& ima, image2d<size_type>& parent,
			const Neighborhood& nbh,
			StrictWeakOrdering cmp, size_type* Send)
	{
	  maxtree_flood_algorithm x(ima, parent, nbh, cmp, Send);
	  if (!parallel)
	    assert((x.m_out + ima.domain().size()) == Send);
	  (void) x;
	}

      private:
	typedef mln_ch_value(I, bool) J;

	StrictWeakOrdering	       m_cmp;
	indexer<V, StrictWeakOrdering> m_h;

	const I&		    m_ima;
	mln_ch_value(I, elt_type)&  m_parent;
	mln_ch_value(I, bool)	    m_deja_vu;
	std::array<typename I::difference_type, Neighborhood::static_size> m_nbh_delta_indexes;

	bounded_hqueue<size_type, nlevels>  m_q;
	bool			m_has_repr[nlevels];
	elt_type		m_repr[nlevels];
	size_type*		m_out;
	size_type		m_first_index;
	size_type		m_last_index;
      };

    } // end of namespace mln::morpho::internal


    // template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
    // image2d<std::size_t>
    // maxtree_hqueue(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
    // {
    //   image2d<std::size_t> parent;
    //   resize(parent, ima);

    //   internal::maxtree_flood_algorithm<image2d<V>, Neighborhood, StrictWeakOrdering>::run(ima, parent, nbh, cmp);
    //   return parent;
    // }

  }  // end of namespace mln::morpho

} // end of namespace mln

#endif // ! MLN_CORE_MORPHO_MAXTREE_HQUEUE_HPP
