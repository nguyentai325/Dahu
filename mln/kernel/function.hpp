#ifndef MLN_KERNEL_FUNCTION_HPP
# define MLN_KERNEL_FUNCTION_HPP

# include <boost/proto/proto.hpp>

namespace mln
{

  namespace kernel
  {

    template <class F>
    using Function = typename boost::proto::terminal<F>::type;

  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_FUNCTION_HPP
