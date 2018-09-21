#ifndef MLN_KERNELV2_DETAILS_CONTEXT_HPP
# define MLN_KERNELV2_DETAILS_CONTEXT_HPP

# include <boost/proto/proto.hpp>
# include <mln/core/image/image.hpp>
# include <mln/kernelv2/details/tags.hpp>
namespace mln
{
  namespace kernel
  {
    namespace details
    {
      namespace proto = boost::proto;

      // The kernel abstract context does nothing
      // but helps to compute the type of expression
      // without needing a real evaluation context
      struct kernel_abstract_context
      {

        template <class Expr, class Tag = typename proto::tag_of<Expr>::type>
        struct eval
          : proto::default_eval<Expr, const kernel_abstract_context, Tag>
        {
        };


        template <class Expr, int k>
        struct eval<Expr, tag::image_call_p<k> >
        {
          typedef typename proto::result_of::value<Expr>::type I;
          typedef mln_reference(I) result_type;

          result_type operator() (Expr, kernel_abstract_context) const;
        };

        template <class Expr, int k>
        struct eval<Expr, tag::image_call_n<k> >
        {
          typedef typename proto::result_of::value<Expr>::type type_;
          typedef typename std::remove_reference<type_>::type::type I;
          typedef mln_reference(I) result_type;

          result_type operator() (Expr, kernel_abstract_context) const;
        };

        template <class Expr>
        struct eval<Expr, tag::aggregate>
        {
          typedef typename proto::result_of::left<Expr>         accu_expr;
          typedef typename proto::result_of::value<accu_expr>   accu;
          typedef typename accu::result_type                    result_type;

          result_type operator() (Expr, kernel_abstract_context) const;
        };
      };

      template <class PVTuple, class NVTuple>
      struct base_context
      {
        PVTuple pvals;
        NVTuple nvals;
      };

    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!MLN_KERNELV2_DETAILS_CONTEXT_HPP
