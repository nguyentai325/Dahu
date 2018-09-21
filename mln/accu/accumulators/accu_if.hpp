#ifndef MLN_ACCU_ACCUMULATORS_ACCU_IF_HPP
# define MLN_ACCU_ACCUMULATORS_ACCU_IF_HPP

# include <mln/accu/accumulator_base.hpp>

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {

      template <class Accu, class Predicate, class ArgType = typename Accu::argument_type>
      struct accu_if :
	accumulator_base< accu_if<Accu, Predicate, ArgType>,
			  ArgType,
			  typename Accu::result_type,
			  typename Accu::feature
			  >
      {
	typedef ArgType			     argument_type;
	typedef typename Accu::result_type   result_type;
	typedef typename Accu::feature	     feature;


	accu_if(const Accu& accu = Accu(),
		const Predicate& pred = Predicate())
	  : m_accu(accu),
	    m_pred(pred)
	{
	}

	void init()
	{
	  m_accu.init();
	}

	void take(const argument_type& arg)
	{
	  if (m_pred(arg))
	    m_accu.take(arg);
	}

	void take(const accu_if& other)
	{
	  m_accu.take(other.m_accu);
	}

	template <class Feat>
	friend
	result_type extract(const accu_if& accu, Feat feat)
	{
	  return extract(accu.m_accu, feat);
	}


      private:
	Accu m_accu;
	Predicate m_pred;
      };


    }

  }

}

#endif // ! ACCU_IF_HPP
