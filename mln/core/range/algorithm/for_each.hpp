#ifndef MLN_CORE_RANGE_ALGORITHM_FOR_EACH_HPP
# define MLN_CORE_RANGE_ALGORITHM_FOR_EACH_HPP

# include <mln/core/range/range.hpp>

namespace mln
{

  namespace range
  {

    template<class InputRange, class Function>
    Function
    for_each(InputRange&& rng,  Function f)
    {
      mln_iter(v, rng);
      mln_forall(v)
	f(*v);
      return f;
    }

  }

}

#endif // ! MLN_CORE_RANGE_ALGORITHM_FOR_EACH_HPP
