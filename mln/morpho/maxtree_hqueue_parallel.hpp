#ifndef MAXTREE_HQUEUE_PARALLEL_HPP
# define MAXTREE_HQUEUE_PARALLEL_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/sub_image.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/wrt_offset.hpp>

# include <mln/io/imprint.hpp>
# include <mln/morpho/maxtree_hqueue.hpp>
# include <mln/morpho/maxtree_routines.hpp>
# include <mln/morpho/merge_tree.hpp>
# include <mln/morpho/canonize.hpp>

# include <tbb/parallel_reduce.h>
# include <tbb/parallel_for.h>
# include <tbb/task_scheduler_init.h>

namespace mln
{

  namespace morpho
  {

    namespace impl
    {

      template <typename V, typename Neighborhood, typename StrictWeakOrdering, bool parallel>
      struct MaxTreeAlgorithmHQ
      {
	typedef typename image2d<V>::size_type size_type;
	static constexpr const size_type UNINITIALIZED = std::numeric_limits<size_type>::max();
	static constexpr const size_type INQUEUE = 0;
	static constexpr const bool use_dejavu = false;

        MaxTreeAlgorithmHQ(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp)
          : m_ima (ima), m_nbh (nbh), m_cmp(cmp), m_has_previous(false)
        {
	  if (!use_dejavu) {
	    resize(m_parent, ima).init((size_type) UNINITIALIZED);
	    extension::fill(m_parent, (size_type) INQUEUE);
	  } else {
	    resize(m_parent, ima);
	  }

	  m_nsplit = 0;
	  if (!parallel) {
	    m_S.resize(ima.domain().size());
	  }
        }


        MaxTreeAlgorithmHQ(MaxTreeAlgorithmHQ& other, tbb::split)
          : m_ima(other.m_ima), m_nbh(other.m_nbh), m_cmp(other.m_cmp), m_parent(other.m_parent),
            m_has_previous(false)
        {
	  m_nsplit = 0;
        }



        void
        operator() (const box2d& domain)
        {
	  image2d<V> ima = m_ima | domain;
          image2d<size_type> parent = m_parent | domain;

	  if (parallel)
	    internal::maxtree_flood_algorithm<image2d<V>, Neighborhood, StrictWeakOrdering, parallel>::run(ima, parent, m_nbh, m_cmp, NULL);
	  else
	    internal::maxtree_flood_algorithm<image2d<V>, Neighborhood, StrictWeakOrdering, parallel>::run(ima, parent, m_nbh, m_cmp, &m_S[0] + domain.size());

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

        void join(MaxTreeAlgorithmHQ& other, bool joindomain = true)
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

        image2d<size_type>	 m_parent;
        bool			 m_has_previous;
        box2d			 m_current_domain;
	unsigned		 m_nsplit;
	std::vector<size_type> m_S;
      };




      namespace parallel
      {
	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_hqueue(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmHQ<V, Neighborhood, StrictWeakOrdering, true> algo(ima, nbh, cmp);
	  int nmaxsplit = tbb::task_scheduler_init::default_num_threads() * 4;
	  int grain = std::max(ima.nrows() / nmaxsplit, 1u);
	  std::cout << "Grain: " << grain << std::endl;
	  tbb::parallel_reduce(grain_box2d(ima.domain(), grain), algo, tbb::auto_partitioner());

	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;
	  image2d<size_type>& parent = algo.m_parent;
	  std::vector<size_type> S(ima.domain().size());
	  canonize(ima, parent, &S[0]);

	  // MaxtreeCanonizationAlgorithm<V> canonizer(ima, parent);
	  // tbb::parallel_for(grain_box2d(ima.domain(), grain), canonizer, tbb::auto_partitioner());

	  return std::make_pair(std::move(parent), std::move(S));
	}
      }


      namespace serial
      {

	template <typename V, typename Neighborhood, typename StrictWeakOrdering = std::less<V> >
	std::pair< image2d<typename image2d<V>::size_type>, std::vector<typename image2d<V>::size_type> >
	maxtree_hqueue(const image2d<V>& ima, const Neighborhood& nbh, StrictWeakOrdering cmp = StrictWeakOrdering())
	{
	  typedef typename image2d<V>::size_type size_type;
	  MaxTreeAlgorithmHQ<V, Neighborhood, StrictWeakOrdering, false> algo(ima, nbh, cmp);
	  algo(ima.domain());
	  std::cout << "Number of split: " << algo.m_nsplit << std::endl;
	  return std::make_pair(std::move(algo.m_parent), std::move(algo.m_S));
	}

      }

    }

  }

}

#endif // !MLN_MORPHO_MAXTREE_HQUEUE_PARALLEL_HPP
