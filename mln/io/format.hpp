#ifndef MLN_IO_FORMAT_HPP
# define MLN_IO_FORMAT_HPP

# include <iostream>
# include <mln/core/grays.hpp>
# include <cmath>
# include <sstream>

namespace mln
{

  namespace io
  {

    template <typename T>
    inline
    std::ostream&
    format(std::ostream& os, const T& v);


    /// \brief Retrieve the maximum width of the a type
    template <class T, class enable = void>
    struct frmt_max_width;

    /****************************/
    /** Implementation         **/
    /****************************/

    template <typename T>
    inline
    std::ostream&
    format(std::ostream& os, const T& v)
    {
      return os << v;
    }

    template <class T, class enable>
    struct frmt_max_width
      : std::integral_constant<int, 0>
    {
      int
      operator() (const T& v) const
      {
	std::ostringstream s;
	format(s, v);
	return s.tellp();
      }

    };

    template <class T>
    struct frmt_max_width<T, typename std::enable_if< std::is_integral<T>::value >::type >
      : std::integral_constant<int, (int) (1 + value_traits<T>::quant / 3.3219280948873622)>
    {

      int
      operator() (const T& v) const
      {
	return (int)(std::log10(std::abs(v))) + 1 ;
      }

    };


    inline
    std::ostream&
    format(std::ostream& os, const uint8& v)
    {
      return os << (int)v;
    }

  }

}

#endif // ! MLN_IO_FORMAT_HPP
