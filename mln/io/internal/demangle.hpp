#ifndef DEMANGLE_HPP
# define DEMANGLE_HPP

# include <string>

# if (__GNUC__ && __cplusplus && __GNUC__ >= 3)
#   define MLN_IO_USE_DEMANGLING
#   include <cxxabi.h>
# endif // __GNUC__

namespace mln
{

  namespace io
  {

    namespace internal
    {


      std::string
      demangle(const char* name);


      /*******************/
      /** Implementation */
      /*******************/


      inline
      std::string
      demangle(const char* name)
      {
        std::string res;

# ifdef MLN_IO_USE_DEMANGLING
        char*       realname;
        int         stat;

        realname = abi::__cxa_demangle(name,NULL,NULL,&stat);
        res = (realname != NULL) ? realname : name;
        free(realname);
# else
        res = name;
# endif
        return res;
      }


    }

  }

}


#endif // ! DEMANGLE_HPP
