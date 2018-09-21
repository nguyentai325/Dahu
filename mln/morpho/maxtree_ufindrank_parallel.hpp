#ifndef MAXTREE_UFINDRANK_PARALLEL_HPP
# define MAXTREE_UFINDRANK_PARALLEL_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/sort_sites.hpp>
# include <mln/core/wrt_offset.hpp>

# include <mln/morpho/merge_tree.hpp>
# include <mln/morpho/canonize.hpp>

# include <tbb/parallel_reduce.h>
# include <tbb/parallel_for.h>


namespace mln
{

  namespace morpho
  {

    namespace impl
    {

      template <typename V, typename Neighborhood, typename StrictWeakOrdering, typename size_type, bool parallel>
      struct MaxTreeAlgorithmUFRank
      {
	static constexpr bool level_compression = false;// value_traits<V>::quant < 18;
	static constexpr bool use_dejavu = false;
	static constexpr size_type UNINITIALIZED = std::numeric_limits<size_type>::max();

	struct aux_data
	{
	  unsigned	rank;
	  size_type	zpar;
	  size_type	repr;
	};


        MaxTreeAlgorithmUFRank(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp)
          : m_ima (ima), m_nbh (nbh), m_cmp(cmp), m_has_previous(false)
        {
          resize(m_parent, ima);
	  aux_data x = {0, UNINITIALIZED, 0};
	  resize(m_aux, ima).init(x);

	  size_type n = m_ima.domain().size();
	  m_S = new std::vector<size_type>(n);
	  m_nsplit = 0;
        }


        MaxTreeAlgorithmUFRank(MaxTreeAlgorithmUFRank& other, tbb::split)
          : m_ima(other.m_ima), m_nbh(other.m_nbh), m_cmp(other.m_cmp), m_parent(other.m_parent),
            m_aux(other.m_aux), m_has_previous(false), m_S (other.m_S)
        {
	  m_nsplit = 0;
        }

	size_type
	zfindroot(image2d<aux_data>& aux, size_type p)
	{
	  size_type q = aux[p].zpar;
	  if (p != q)
	    return aux[p].zpar = zfindroot(aux, q);
	  else
	    return q;
	}


	void
	unionfind(const grain_box2d& domain)
	{
          image2d<V> ima = m_ima | domain;
          image2d<size_type>  parent = m_parent | domain;
          image2d<aux_data>  aux = m_aux | domain;
          image2d<bool> deja_vu;

	  if (use_dejavu)
	    resize(deja_vu, ima).init(false);

	  const box2d& d = m_ima.domain();
	  size_type i = (d.pmax[1] - d.pmin[1]) * (domain.pmin[0] - d.pmin[0]);
	  size_type* S = m_S->data() + i;

	  //int djvu_offset = deja_vu.index_of_point(domain.pmin) - ima.index_of_point(domain.pmin);
          sort_indexes_it(ima, S, m_cmp);
	  auto dindexes = wrt_delta_index(ima, m_nbh.dpoints);
	  size_type* R = S + ima.domain().size();

	  size_type first_index = ima.index_of_point(ima.domain().pmin);
	  size_type last_index = ima.index_of_point(ima.domain().pmax) - ima.index_strides()[0];

          for (int i = ima.domain().size()-1; i >= 0; --i)
            {
	      size_type p = S[i];
              // make set
              {
                parent[p] = p;
                aux[p].zpar = p;
		aux[p].repr = p;
		if (use_dejavu)
		  deja_vu[p] = true;
              }

	      size_type z = p; // zpar of p
	      size_type x = p; // aux[z].repr
              for (unsigned j = 0; j < dindexes.size(); ++j)
		{
		  size_type n = p + dindexes[j];

		  bool processed;
		  if (use_dejavu)
		    processed = deja_vu[n];
		  else if (!parallel) // no bound checking
		    processed = (aux[n].zpar != UNINITIALIZED);
		  else
		    processed = (aux[n].zpar != UNINITIALIZED) and (first_index <= n and n < last_index);

		  if (processed)
                  {
		    size_type r = zfindroot(aux, n);
                    if (r != z) // make union
                      {
			if (level_compression and ima[aux[r].repr] == ima[x]) // FIXME: use !(b>a) for a==b because ima[x]<=ima[repr[r]]
			  std::swap(x, aux[r].repr);
			parent[aux[r].repr] = x;

			if (level_compression and !parallel)
			  *(--R) = aux[r].repr;

			if (aux[z].rank < aux[r].rank) { //we merge z to r
			  aux[z].zpar = r;
			  aux[r].repr = x;
			  z = r;
			} else if (aux[z].rank > aux[r].rank) { // merge r to z
			  aux[r].zpar = z;
			  aux[z].repr = x;
			} else { // same height
			  aux[r].zpar = z;
			  aux[z].repr = x;
			  aux[z].rank += 1;
			}
                      }
                  }
		}
            }
	  if (level_compression and !parallel) {
	    *(--R) = parent[S[0]];
	    assert(R == S);
	    check_S(parent, S, S + ima.domain().size());
	  }
	}

        void
        operator() (const box2d& domain)
        {
	  //std::cout << domain << std::endl;
	  unionfind(domain);

          if (m_has_previous)
	    {
	      this->join(*this, false);
	      m_current_domain.join(domain);
	      m_nsplit += 1;
	    }
	  else
	    {
	      m_current_domain = domain;
	      m_has_previous = true;
	    }
        }

        void join(MaxTreeAlgorithmUFRank& other, bool joindomain = true)
        {
	  mln_precondition(m_has_previous);

          merge_tree(m_ima, m_parent, this->m_current_domain, m_cmp);

	  if (joindomain)
	    {
	      m_current_domain.join(other.m_current_domain);
	      m_nsplit += other.m_nsplit + 1;
	    }
        }


      public:
        const image2d<V>&  m_ima;
        Neighborhood	   m_nbh;
        StrictWeakOrdering m_cmp;

        image2d<size_type> m_parent;
        image2d<aux_data>	     m_aux;

        bool	          m_has_previous;
        box2d	          m_current_domain;

        std::vector<size_type>* m_S;
	unsigned	     m_nsplit;
      };


      namespace parallel
      {
	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair<image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_ufindrank(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmUFRank<V, Neighborhood, StrictWeakOrdering, size_type, true> algo(ima, nbh, cmp);
	  int grain = std::max(ima.nrows() / 64, 1u);
	  std::cout << "Grain: " << grain << std::endl;
	  tbb::parallel_reduce(grain_box2d(ima.domain(), grain), algo, tbb::auto_partitioner());

	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;

	  image2d<size_type>& parent = algo.m_parent;
	  std::vector<size_type> S = std::move(*(algo.m_S));
	  delete algo.m_S;

	  canonize(ima, parent, &S[0]);
	  return std::make_pair(std::move(parent), std::move(S));
	}

      }


      namespace serial
      {

	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair<image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_ufindrank(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmUFRank<V, Neighborhood, StrictWeakOrdering, size_type, false> algo(ima, nbh, cmp);
	  algo(ima.domain());

	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;

	  std::vector<size_type> S = std::move(*(algo.m_S));
	  delete algo.m_S;

	  // canonization
	  // image2d<std::size_t>& parent = algo.m_parent;
	  // for (std::size_t p : S)
	  //   {
	  //     std::size_t q = parent[p];
	  //     if (ima[parent[q]] == ima[q])
	  // 	parent[p] = parent[q];
	  //   }
	  canonize(ima, S, algo.m_parent);
	  return std::make_pair(std::move(algo.m_parent), std::move(S));
	}

      }

    }


    // template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
    // image2d<point2d>
    // maxtree_ufind_parallel(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
    // {
    //   impl::MaxTreeAlgorithmUF<V, Neighborhood, StrictWeakOrdering> algo(ima, nbh, cmp);
    //   tbb::parallel_reduce(ima.domain(), algo, tbb::auto_partitioner());

    //   // canonization
    //   image2d<point2d>& parent = algo.m_parent;
    //   mln_foreach(auto& p, parent.values())
    // 	p = internal::zfind_repr(ima, parent, p);

    //   return algo.m_parent;
    // }



  }

}

#endif // !MLN_MORPHO_MAXTREE_UFINDRANK_PARALLEL_HPP
