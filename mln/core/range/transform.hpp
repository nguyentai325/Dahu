#ifndef MLN_CORE_RANGE_TRANSFORM_HPP
# define MLN_CORE_RANGE_TRANSFORM_HPP

# include <mln/core/iterator/transform_iterator.hpp>

namespace mln
{

  template <typename InputRange, typename UnaryFunction>
  struct transformed_range
  {
  private:
    template <typename R, typename F>
    friend struct transformed_range;

    typedef typename std::remove_reference<InputRange>::type R;

  public:
    typedef transform_iterator<typename R::iterator, UnaryFunction>		iterator;
    typedef transform_iterator<typename R::const_iterator, UnaryFunction>	const_iterator;
    typedef typename iterator::value_type						value_type;
    typedef typename iterator::reference						reference;

    transformed_range(InputRange&& rng, const UnaryFunction& fun)
      : m_rng (std::forward<InputRange>(rng)), m_fun(fun)
    {
    }

    transformed_range(const transformed_range&) = default;

    //interop
    template <typename R2, typename F2>
    transformed_range(const transformed_range<R2, F2>& other,
		      typename std::enable_if<std::is_convertible<R2, InputRange>::value &&
		      std::is_convertible<F2, UnaryFunction>::value>::type* = NULL)
      : m_rng(other.m_rng),
	m_fun(other.m_fun)
    {
    }

    const_iterator iter() const
    {
      return const_iterator( m_rng.iter(), m_fun );
    }

    iterator iter()
    {
      return iterator( m_rng.iter(), m_fun );
    }

    std::size_t size() const
    {
      return m_rng.size();
    }

  private:
    InputRange	m_rng;
    UnaryFunction     m_fun;
  };



  namespace rng
  {

    template <typename InputRange, typename UnaryFunction>
    transformed_range<InputRange, UnaryFunction>
    transform(InputRange&& rng, const UnaryFunction& fun)
    {
      return transformed_range<InputRange, UnaryFunction>(rng, fun);
    }

  } // end of namespace mln::rng

} // end of namespace mln

#endif // ! MLN_CORE_RANGE_TRANSFORM_HPP
