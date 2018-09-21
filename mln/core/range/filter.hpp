#ifndef MLN_CORE_RANGE_FILTER_HPP
# define MLN_CORE_RANGE_FILTER_HPP

# include <mln/core/iterator/filter_iterator.hpp>

namespace mln
{

    template <typename InputRange, typename Predicate>
    struct filtered_range
    {
    private:
      typedef typename std::remove_reference<InputRange>::type R;

    public:
      typedef filter_iterator<typename R::iterator, Predicate>		iterator;
      typedef filter_iterator<typename R::const_iterator, Predicate>	const_iterator;
      typedef typename iterator::value_type				value_type;
      typedef typename iterator::reference				reference;

      filtered_range(InputRange&& rng, const Predicate& pred)
	: m_rng (std::forward<InputRange>(rng)), m_pred(pred)
      {
      }

      const_iterator iter() const
      {
	return const_iterator( m_rng.iter(), m_pred );
      }

      iterator iter()
      {
	return iterator( m_rng.iter(), m_pred );
      }

      bool has(const value_type& v) const
      {
	return m_pred(v) and m_rng.has(v);
      }

    private:
      InputRange  m_rng;
      Predicate   m_pred;
    };



  namespace rng
  {

    template <typename InputRange, typename Predicate>
    filtered_range<InputRange, Predicate>
    filter(InputRange&& rng, Predicate pred)
    {
      return filtered_range<InputRange, Predicate>(std::forward<InputRange>(rng), pred);
    }

  } // end of namespace mln::rng

} // end of namespace mln

#endif // ! MLN_CORE_RANGE_FILTER_HPP
