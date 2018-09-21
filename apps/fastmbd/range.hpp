#ifndef RANGE_HPP
#define RANGE_HPP


# include <mln/io/format.hpp>


namespace mln
{


    template <class V>
    struct range
    {
      range() = default;
      range(const V& v) : lower (v), upper(v) {}
      range(const V& lower_, const V& upper_) : lower (lower_), upper(upper_) {}

      V lower, upper;
    };


    template <typename V>
    std::ostream&
    operator <<(std::ostream& os, const range<V>& rng)
    {
      using namespace mln::io;
      os << '['; format(os, rng.lower);
      os << ','; format(os, rng.upper);
      return os << ']';
    }

}

#endif // RANGE_HPP
