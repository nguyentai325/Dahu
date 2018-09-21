#ifndef MLN_KERNELV2_KERNEL_HPP
# define MLN_KERNELV2_KERNEL_HPP

# include <mln/core/config.hpp>
# include <mln/kernelv2/types.hpp>
# include <mln/kernelv2/execute.hpp>
# include <mln/kernelv2/execute_incremental.hpp>
# include <mln/kernelv2/details/deep_copy_fix.hpp>

namespace mln
{

  namespace kernel
  {

    /// \brief Declare a kernel expression
    template <class Expr>
    typename proto::result_of::deep_copy<Expr>::type
    declare(Expr&& e);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/


    template <class Expr>
    typename proto::result_of::deep_copy<Expr>::type
    declare(Expr&& e)
    {
      return proto::deep_copy(std::forward<Expr>(e));
    }


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNELV2_KERNEL_HPP
