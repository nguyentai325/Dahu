#ifndef MLN_KERNELV2_DETAILS_INNER_CONTEXT_HPP
# define MLN_KERNELV2_DETAILS_INNER_CONTEXT_HPP

# include <mln/kernelv2/details/eval_context.hpp>

namespace mln
{
  namespace kernel
  {
    namespace details
    {

      /// \brief a proto transform that calls init on every
      /// accumulator
      template <class F, class PVal, class NVal>
      struct inner_context;

      struct inner_context_init;
      struct accu_take;
      struct accu_untake;

      template <class PVal, class NVal>
      using inner_context_take = inner_context<accu_take, PVal, NVal>;

      template <class PVal, class NVal>
      using inner_context_untake = inner_context<accu_untake, PVal, NVal>;


      /******************************************/
      /****          Implementation          ****/
      /******************************************/

      struct accu_take
      {
        template <class Accu, class Varg>
        void
        operator() (Accu& x, Varg&& v) const
        {
          x.take(std::forward<Varg>(v));
        }
      };

      struct accu_untake
      {
        template <class Accu, class Varg>
        void
        operator() (Accu& x, Varg&& v) const
        {
          x.untake(std::forward<Varg>(v));
        }
      };


      template <class F, class PVal, class NVal>
      struct inner_context
      {
        template <typename Expr, typename Tag = typename proto::tag_of<Expr>::type>
        struct eval : proto::null_eval<Expr, inner_context<F, PVal, NVal> >
        {
        };

        template <typename Expr>
        struct eval<Expr, tag::aggregate>
        {
          typedef void result_type;

          void
          operator () (Expr& x, inner_context& ctx) const
          {
            auto& accu = proto::value(proto::left(x));
            eval_context<PVal, NVal> ctx2 = {ctx.pvals, ctx.nvals};

            F() (accu, proto::eval(proto::right(x), ctx2));
          }
        };

        PVal pvals;
        NVal nvals;
      };


      struct inner_context_init
      {
        template <typename Expr, typename Tag = typename proto::tag_of<Expr>::type>
        struct eval : proto::null_eval<Expr, inner_context_init>
        {
        };

        template <typename Expr>
        struct eval<Expr, tag::aggregate>
        {
          typedef void result_type;

          void
          operator () (Expr& x, inner_context_init) const
          {
            proto::value(proto::left(x)).init();
          }
        };
      };

    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!MLN_KERNELV2_DETAILS_INNER_CONTEXT_HPP
