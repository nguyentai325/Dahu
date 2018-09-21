#ifndef MLN_CORE_RANGE_ALGORITHM_GENERATE_HPP
# define MLN_CORE_RANGE_ALGORITHM_GENERATE_HPP

# include <mln/core/range/range.hpp>

namespace mln
{
  namespace range
  {

    template <class Range, class Generator>
    void generate(Range&& rng, Generator f)
    {
      mln_iter(vin, rng);
      mln_forall(vin)
        *vin = f();
    }

  } // end of namespace mln::range

} // end of namespace mln

#endif //!MLN_CORE_RANGE_ALGORITHM_GENERATE_HPP
