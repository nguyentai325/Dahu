#ifndef MAXTREE_UFIND_PARALLEL_HPP
# define MAXTREE_UFIND_PARALLEL_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/sort_indexes.hpp>
# include <mln/core/algorithm/fill.hpp>
# include <mln/core/wrt_offset.hpp>

# include <mln/io/imprint.hpp>
# include <mln/morpho/merge_tree.hpp>
# include <mln/morpho/canonize.hpp>
# include <alloca.h>

# include <tbb/parallel_reduce.h>
# include <tbb/parallel_for.h>

# include <tbb/mutex.h>

namespace mln
{

  namespace morpho
  {

    namespace impl
    {



      template <typename V, typename Neighborhood, typename StrictWeakOrdering, typename size_type, bool parallel>
      struct MaxTreeAlgorithmUF
      {
#  ifdef LEVEL_COMPRESSION
	static constexpr bool level_compression =  true;
#  elif defined(NO_LEVEL_COMPRESSION)
	static constexpr bool level_compression =  false;
#  else
	static constexpr bool level_compression = value_traits<V>::quant < 18;
#  endif

	static constexpr bool use_dejavu = false;
	static constexpr size_type UNINITIALIZED = std::numeric_limits<size_type>::max();

        MaxTreeAlgorithmUF(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp)
          : m_ima (ima), m_nbh (nbh), m_cmp(cmp), m_has_previous(false)
        {
	  size_type n = m_ima.domain().size();
	  m_S = new std::vector<size_type>(n);

          resize(m_parent, ima);
          resize(m_zpar, ima).init(UNINITIALIZED);
	  m_nsplit = 0;
        }


        MaxTreeAlgorithmUF(MaxTreeAlgorithmUF& other, tbb::split)
          : m_ima(other.m_ima), m_nbh(other.m_nbh), m_cmp(other.m_cmp), m_parent(other.m_parent),
            m_zpar(other.m_zpar), m_has_previous(false), m_S (other.m_S)
        {
	  m_nsplit = 0;
        }

	static
	inline
	size_type
	zfindroot(image2d<size_type>& parent, const size_type& p)
	{
	  if (parent[p] == p)
	    return p;
	  else
	    return parent[p] = zfindroot(parent, parent[p]);
	}

	void
	unionfind(const box2d& domain)
	{
          image2d<V> ima = m_ima | domain;
          image2d<size_type> parent = m_parent | domain;
          image2d<size_type> zpar = m_zpar | domain;


	  const box2d& d = m_ima.domain();
	  size_type i = (d.pmax[1] - d.pmin[1]) * (domain.pmin[0] - d.pmin[0]);
	  size_type* S = m_S->data() + i;

	  image2d<bool> deja_vu;
	  if (use_dejavu)
	    resize(deja_vu, ima).init(false);

	  size_type first_index = ima.index_of_point(ima.domain().pmin);
	  size_type last_index = ima.index_of_point(ima.domain().pmax) - ima.index_strides()[0];


          sort_indexes_it(ima, S, m_cmp);
	  auto dindexes = wrt_delta_index(ima, m_nbh.dpoints);


          size_type* R = S + ima.domain().size();
          for (int i = ima.domain().size() - 1; i >= 0; --i)
            {
	      size_type p = S[i];
	      //std::cout << deja_vu.point_at_index(djvu_offset + p) << std::endl;
              // make set
	      assert(!use_dejavu or domain.has(deja_vu.point_at_index(p)));
              {
                parent[p] = p;
                zpar[p] = p;
		if (use_dejavu)
		  deja_vu[p] = true;
              }

	      size_type z = p; // zpar of p
              for (unsigned j = 0; j < dindexes.size(); ++j)
		{
		  size_type n = p + dindexes[j];

		  bool processed;
		  if (use_dejavu)
		    processed = deja_vu[n];
		  else if (!parallel) // no bound checking
		    processed = (zpar[n] != UNINITIALIZED);
		  else
		    processed = (zpar[n] != UNINITIALIZED) and (first_index <= n and n < last_index);

		  if (processed)
                  {
		    size_type r = zfindroot(zpar, n);

                    if (r != z) // make union
                      {
			if (level_compression and ima[r] == ima[z])
			  std::swap(r, z);

			parent[r] = z;
			zpar[r] = z;
			if (level_compression and !parallel)
			  *(--R) = r;
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
	unionfind_line(int row)
	{
	  static constexpr std::size_t nvalues = 1ul << value_traits<V>::quant;
	  point2d p_ = m_ima.domain().pmin;
	  p_[0] = row;

	  int ncols = m_ima.ncols();
	  int sz = 0;
	  size_type prec = m_ima.index_of_point(p_);
	  size_type p = prec + 1;
	  size_type* stack;

	  if (value_traits<V>::quant <= 16)
	    stack = (size_type*) alloca(nvalues * sizeof(size_t));
	  else
	    stack = new size_type[nvalues];

	  for (int i = 1; i < ncols; ++i, ++p)
	    {
	      if (m_cmp(m_ima[prec],m_ima[p])) // m_ima[prec] < m_ima[p] => start new component
		{
		  stack[sz++] = prec;
		  //m_parent[prec] = prec;
		  prec = p;
		}
	      else if (not m_cmp(m_ima[p], m_ima[prec])) // m_ima[p] == m_ima[prec] => extend component
		{
		  m_parent[p] = prec;
		}
	      else // m_ima[p] < m_ima[prec] => we need to attach prec to its m_parent
		{
		  while (sz > 0 and not m_cmp(m_ima[stack[sz-1]], m_ima[p]))
		    {
		      m_parent[prec] = stack[sz-1];
		      prec = stack[sz-1];
		      --sz;
		    }
		  // we have m_ima[p] <= m_ima[prec]
		  if (m_cmp(m_ima[p], m_ima[prec])) // m_ima[p] < m_ima[prec] => attach prec to p, p new component
		    {
		      m_parent[prec] = p;
		      prec = p;
		    }
		  else                        // m_ima[p] == m_ima[prec] => attach p to prec (canonization)
		    {
		      m_parent[p] = prec;
		    }
		}
	    }

	  // Attach last point (i.e prec)
	  while (sz > 0)
	    {
	      m_parent[prec] = stack[sz-1];
	      prec = stack[sz-1];
	      --sz;
	    }
	  m_parent[prec] = prec;

	  if (value_traits<V>::quant > 16)
	    delete [] stack;
	}


        void
        operator() (const box2d& domain)
        {
	  //std::cout << domain << std::endl;
	  if (domain.shape()[0] > 1)
	    unionfind(domain);
	  else
	    unionfind_line(domain.pmin[0]);


          if (m_has_previous)
	    {
	      merge_tree(m_ima, m_parent, this->m_current_domain, m_cmp);
	      m_current_domain.join(domain);

	      // const box2d& d = m_ima.domain();
	      // size_type w = (d.pmax[1] - d.pmin[1]);
	      // size_type begin1 = w * (m_current_domain.pmin[0] - d.pmin[0]);
	      // size_type end1 = w * (m_current_domain.pmax[0] - d.pmin[0]);
	      // size_type end2 = w * (domain.pmax[0] - d.pmin[0]);
              // check_S(m_parent, m_Ssrc->data() + begin1, m_Ssrc->data() + end1);
              // check_S(m_parent, m_Ssrc->data() + end1, m_Ssrc->data() + end2);

	      // fill(m_dejavu | m_current_domain, false);
	      // merge_S(m_parent, m_dejavu, m_Ssrc->data() + begin1, m_Ssrc->data() + end1,
	      //         m_Ssrc->data() + end1, m_Ssrc->data() + end2, m_Sdst->data() + begin1);
	      // std::swap(*m_Ssrc, *m_Sdst);
	      m_nsplit += 1;
	    }
	  else
	    {
	      m_current_domain = domain;
	      m_has_previous = true;
	    }
        }

        void join(MaxTreeAlgorithmUF& other)
        {
	  mln_precondition(m_has_previous);

	  // Merge trees
          merge_tree(m_ima, m_parent, this->m_current_domain, m_cmp);
          m_current_domain.join(other.m_current_domain);
          m_nsplit += other.m_nsplit + 1;
        }


      public:
        const image2d<V>&  m_ima;
        Neighborhood	   m_nbh;
        StrictWeakOrdering m_cmp;

        image2d<size_type> m_parent;
        image2d<size_type> m_zpar;

        bool	             m_has_previous;
        box2d	             m_current_domain;

	std::vector<size_type>*  m_S;

	unsigned	     m_nsplit;
      };

      template <typename V, typename Neighborhood, typename StrictWeakOrdering, typename size_type, bool parallel>
      constexpr size_type MaxTreeAlgorithmUF<V, Neighborhood, StrictWeakOrdering, size_type, parallel>::UNINITIALIZED;


      namespace parallel
      {
	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_ufind(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmUF<V, Neighborhood, StrictWeakOrdering, size_type, true> algo(ima, nbh, cmp);
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

	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_ufind_line(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmUF<V, Neighborhood, StrictWeakOrdering, size_type, true> algo(ima, nbh, cmp);
	  std::cout << "Grain: " <<  1 << std::endl;
	  tbb::parallel_reduce(ima.domain(), algo, tbb::simple_partitioner());

	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;


	  image2d<size_type>&    parent = algo.m_parent;
	  std::vector<size_type> S = std::move(*(algo.m_S));
	  delete algo.m_S;

	  canonize(ima, parent, &S[0]);
	  return std::make_pair(std::move(parent), std::move(S));
	}


      }


      namespace serial
      {

	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_ufind(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmUF<V, Neighborhood, StrictWeakOrdering, size_type, false> algo(ima, nbh, cmp);
	  algo(ima.domain());

	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;
	  // canonization
	  // image2d<size_type>& parent = algo.m_parent;
	  // for (size_type p : algo.m_S)
	  //   {
	  //     size_type q = parent[p];
	  //     if (ima[parent[q]] == ima[q])
	  // 	parent[p] = parent[q];
	  //   }

	  std::vector<size_type> S = std::move(*(algo.m_S));
	  delete algo.m_S;

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

#endif // !MLN_MORPHO_MAXTREE_UFIND_PARALLEL_HPP
