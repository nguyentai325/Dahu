#ifndef ACCU_CVXHULL_HPP
# define ACCU_CVXHULL_HPP

# include <mln/accu/accumulator.hpp>
# include <mln/accu/accumulator_base.hpp>
# include <mln/accu/accumulators/cvxhull_impl.hpp>
namespace mln
{

  namespace accu
  {
    //FWD
    namespace accumulators
    {
      template <typename P>
      struct cvxhull;
    }


    namespace features
    {

      struct cvxhull : simple_feature<cvxhull>
      {
	template <typename P>
	struct apply
	{
	  typedef accumulators::cvxhull<P> type;
	};
      };

      struct pvector : simple_feature<pvector>
      {
	template <typename P>
	struct apply
	{
	  typedef accumulators::cvxhull<P> type;
	};
      };

    }

    namespace extractor
    {
      template <typename A>
      auto
      cvxhull(const Accumulator<A>& accu)
	-> decltype( extract(accu, features::cvxhull ()) )
      {
	return extract(accu, features::cvxhull ());
      }

      template <typename A>
      auto
      pvector(const Accumulator<A>& accu)
	-> decltype( extract(accu, features::pvector ()) )
      {
	return extract(accu, features::pvector ());
      }

    }


    namespace accumulators
    {

      template <typename P>
      struct cvxhull : accumulator_base< cvxhull<P>, P, std::vector<P>, features::cvxhull >
      {
	typedef P				argument_type;
	typedef std::vector<P>			result_type;
	typedef features::cvxhull		feature;
	typedef boost::mpl::set<feature>	provides;

	cvxhull()
	{
	}


	void init()
	{
	  m_points.clear();
	}

	void take(const P& p)
	{
	  mln_precondition(m_points.empty() || m_points.back() < p);
	  m_points.push_back(p);
	}

	template <typename A>
	void take(const Accumulator<A>& other)
	{
	  std::vector<P>& ovec = extractor::pvector(other);
	  std::vector<P> aux(m_points.size() + ovec.size());
	  std::merge(m_points.begin(), m_points.end(), ovec.begin(), ovec.end(), aux.begin());
	  m_points = aux;
	}

	friend
	std::vector<P>&
	extract(cvxhull& accu, features::pvector)
	{
	  return accu.m_points;
	}

	friend
	std::vector<P>
	extract(const cvxhull& accu, features::cvxhull)
	{
	  std::vector<P> pts = convexhull(accu.m_points);
	  return pts;
	}

      private:
	std::vector<P> m_points;
      };

    }

  }

}

#endif // ! ACCU_CVXHULL_HPP
