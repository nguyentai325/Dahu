#undef mln_assertion
#undef mln_precondition
#undef mln_postcondition
#undef MLN_HAS_DEBUG

#include <cassert>

#ifndef _MLN_HAS_NDEBUG_FIRST_
# ifdef NDEBUG
#  define MLN_NDEBUG
# endif
# define _MLN_HAS_NDEBUG_FIRST_
#endif



#ifndef MLN_NDEBUG
# define mln_assertion(expr) assert(expr)
# define MLN_HAS_DEBUG 1
# define MLN_EVAL_IF_DEBUG(expr) expr
#else
# define mln_assertion(expr) (void (0))
# define MLN_HAS_DEBUG 0
# define MLN_EVAL_IF_DEBUG(expr)
#endif

#define mln_precondition(expr) mln_assertion(expr)
#define mln_postcondition(expr) mln_assertion(expr)
