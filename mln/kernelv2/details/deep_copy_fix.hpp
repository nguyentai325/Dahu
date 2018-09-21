#include <boost/proto/deep_copy.hpp>

namespace boost
{
  namespace proto
  {
    namespace detail
    {

    template<class I, int k>
    struct deep_copy_impl< mln::kernel::details::image_call_p_expr<I, k>, 0>
    {
      typedef mln::kernel::details::image_call_n_expr<I,k> Expr;

      typedef
      typename base_expr<
        typename Expr::proto_domain
        , mln::kernel::tag::image_call_p<k>
        , term<typename term_traits<typename Expr::proto_child0>::value_type>
        >::type
      expr_type;

      typedef typename Expr::proto_generator proto_generator;
      typedef typename proto_generator::template result<proto_generator(expr_type)>::type result_type;

      template<typename Expr2, typename S, typename D>
      result_type operator()(Expr2 const &e, S const &, D const &) const
      {
        return proto_generator()(expr_type::make(e.proto_base().child0));
      }
    };

    template<class I, int k>
    struct deep_copy_impl< mln::kernel::details::image_call_n_expr<I,k>, 0>
    {
      typedef mln::kernel::details::image_call_n_expr<I,k> Expr;

      typedef
      typename base_expr<
        typename Expr::proto_domain
        , mln::kernel::tag::image_call_n<k>
        , term<typename term_traits<typename Expr::proto_child0>::value_type>
        >::type
      expr_type;

      typedef typename Expr::proto_generator proto_generator;
      typedef typename proto_generator::template result<proto_generator(expr_type)>::type result_type;

      template<typename Expr2, typename S, typename D>
      result_type operator()(Expr2 const &e, S const &, D const &) const
      {
        return proto_generator()(expr_type::make(e.proto_base().child0));
      }
    };


    }

  }

}
