#ifndef MLN_KERNEL_AGGREGATE_INF_HPP
# define MLN_KERNEL_AGGREGATE_INF_HPP

# include <mln/kernelv2/types.hpp>
# include <mln/accu/accumulators/infsup.hpp>

namespace mln
{
  namespace kernel
  {
    namespace aggregate
    {

      template <class Compare = void>
      using Inf_t = Aggregate< accu::features::inf<Compare> >;


      template <class Compare = void, class Expr>
      auto
      Inf(Expr&& expr)
        -> decltype( Inf_t<Compare> () (std::declval<Expr> () ) )
      {
        return Inf_t<Compare>() (std::forward<Expr>(expr));
      }


    } // end of namespace mln::kernel::aggregate

  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_AGGREGATE_INF_HPP
