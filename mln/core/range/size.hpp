#ifndef RANGE_SIZE_HPP
# define RANGE_SIZE_HPP


# include <type_traits>
# include <mln/core/range/range.hpp>

namespace mln
{

  namespace rng
  {

    template <typename Range>
    std::size_t size(const Range& rng);

    /***********************/
    /*** Specialization  ***/
    /***********************/

    template <typename T, std::size_t N>
    struct array;


    template <typename T, std::size_t N>
    constexpr
    std::size_t size(const std::array<T, N>&)
    {
      return N;
    }

    /***********************/
    /*** Implementation  ***/
    /***********************/

    namespace internal
    {

      template <typename R, typename = void>
      struct size_lookup
      {
	static
	std::size_t size(const R& rng)
	{
	  std::size_t sz = 0;
	  auto it = rng::iter(rng);
	  for (it.init(); !it.finished(); it.next())
	    ++sz;
	  return sz;
	}
      };


      template <typename R>
      struct size_lookup<R, typename std::enable_if<true, decltype(std::declval<R>().size(), (void)0)>::type>
      {
	static
	std::size_t size(const R& rng)
	{
	  return rng.size();
	}
      };

    }

    template <typename Range>
    inline
    std::size_t size(const Range& rng)
    {
      return internal::size_lookup<Range>::size(rng);
    }


  }

}

#endif // ! RANGE_SIZE_HPP
