#ifndef MLN_IO_EXCEPTION_HPP
# define MLN_IO_EXCEPTION_HPP

# include <exception>

namespace mln
{

  namespace io
  {

    class MLNIOException : public std::runtime_error
    {
    public:
      MLNIOException(const std::string& what_arg) :
	std::runtime_error(what_arg)
      {
      }
    };

  }

}
#endif
