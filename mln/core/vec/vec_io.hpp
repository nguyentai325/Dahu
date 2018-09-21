#ifndef MLN_CORE_VEC_VEC_IO_HPP
# define MLN_CORE_VEC_VEC_IO_HPP

# include <mln/core/vec.hpp>
# include <mln/io/format.hpp>
# include <sstream>

namespace mln
{

  namespace internal
  {

    template <typename T, unsigned dim, typename tag>
    inline
    std::ostream&
    operator<< (std::ostream& os, const vec_base<T, dim, tag>& x)
    {
      std::ostringstream s;
      s << '[';
      for (unsigned i = 0; i < dim-1; ++i)
	s << x[i] << ',';
      s << x[dim-1] << ']';

      os << s.str();
      return os;
    }


    template <typename T, unsigned dim, typename tag>
    inline
    std::ostream&
    format (std::ostream& os, const vec_base<T, dim, tag>& x)
    {
      using io::format;
      std::ostringstream s;

      s << '[';
      for (unsigned i = 0; i < dim-1; ++i) {
	format(s, x[i]); s << ',';
      }
      format(s, x[dim-1]); s << ']';

      os << s.str();
      return os;
    }

  }

}

#endif // ! MLN_CORE_VEC_VEC_IO_HPP
