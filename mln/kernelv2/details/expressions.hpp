#ifndef MLN_KERNEL_DETAILS_EXPRESSIONS_HPP
# define MLN_KERNEL_DETAILS_EXPRESSIONS_HPP

# include <boost/proto/proto.hpp>
# include <mln/kernelv2/details/context.hpp>
# include <mln/accu/accumulator.hpp>

namespace mln
{
  namespace kernel
  {
    namespace details
    {
      namespace proto = boost::proto;

      /**
       ** We extends proto with few expression in the DSL
       **
       **
       */

      struct dummy_t {};

      typedef proto::expr<tag::point, proto::term<dummy_t> >     point_expr;
      typedef proto::expr<tag::neighbor, proto::term<dummy_t> >  neighbor_expr;

      template <class I, int k>
      using image_call_p_expr = proto::expr<tag::image_call_p<k>,
                                            proto::term<I> >;

      template <class I, int k>
      using image_call_n_expr = proto::expr<tag::image_call_n<k>,
                                            proto::term<I> >;



      template <class Feature, class Expr>
      using aggregate_expr =  typename proto::result_of::make_expr<
        kernel::tag::aggregate,
        typename mln::accu::accu_of<
          Feature,
          typename std::decay<
            typename proto::result_of::eval<Expr,kernel_abstract_context>::type
            >::type
          >::type,
        Expr>::type;

    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!_MLN_KERNEL_DETAILS_EXPRESSIONS_HPP
