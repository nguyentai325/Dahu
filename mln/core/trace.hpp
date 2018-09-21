#ifndef TRACE_HPP
# define TRACE_HPP

//# include <boost/current_function.hpp>
# include <ctime>
# include <cstdlib>
# include <iostream>
# include <stack>
# include <string>
//# include <mutex>
# include <thread>
# include <tbb/compat/thread>
# include <tbb/combinable.h>
# include <tbb/tick_count.h>

# define mln_entering(NAME) \
  mln::trace::entering(NAME);

# define mln_exiting()  mln::trace::exiting();

namespace mln
{

  namespace trace
  {
    struct trace_t
    {
      explicit trace_t(const std::string& fun)
	: fname (fun)
      {
	clock = tbb::tick_count::now();
      }

      std::string       fname;
      tbb::tick_count   clock;
    };

    static tbb::combinable< std::stack<trace_t> > callstacks;
    //static std::stack<trace_t> callstack;
    static bool verbose = (std::getenv("TRACE") != NULL);
    //static std::mutex stack_mutex;

    static inline
    void entering(const std::string& fname)
    {
      if (verbose) {
        //std::lock_guard<std::mutex> lock(stack_mutex);
        std::stack<trace_t>& callstack = callstacks.local();
        callstack.emplace(fname);

	std::clog << std::string(callstack.size(), ' ')
                  << "#" << std::this_thread::get_id() << " - "
                  << fname << std::endl;
      }
    };


    static inline
    void exiting()
    {
      if (verbose) {
        //std::lock_guard<std::mutex> lock(stack_mutex);
        std::stack<trace_t>& callstack = callstacks.local();
	trace_t tr = callstack.top();
	std::clog << std::string(callstack.size(), ' ')
                  << "#" << std::this_thread::get_id() << " - "
                  << tr.fname
		  << " in " << (tbb::tick_count::now() - tr.clock).seconds() << std::endl;

        callstack.pop();
      }
    };

    inline
    void warn(const std::string& msg)
    {
      if (verbose)
        std::clog << std::string(callstacks.local().size(), ' ')
                  << "#" << std::this_thread::get_id() << " - "
                  << msg << std::endl;
    }

  }

}

#endif // ! TRACE_HPP
