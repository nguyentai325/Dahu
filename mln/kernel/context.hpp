#ifndef MLN_KERNEL_CONTEXT_HPP
# define MLN_KERNEL_CONTEXT_HPP

# include <tuple>
# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/accu/accumulator.hpp>
# include <mln/core/internal/intseq.hpp>
# include <boost/proto/proto.hpp>
# include <mln/kernel/types.hpp>
# include <mln/kernel/aggregate.hpp>
# include <mln/kernel/intro.hpp>
# include <mln/kernel/details/transform.hpp>

namespace mln
{

  namespace kernel
  {

    namespace proto = boost::proto;

    ///
    ///
    ///
    template <class V, class V2,
              class SubexprTypeList,
              class AccuList>
    struct kernel_context;


    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace internal
    {
      struct nil {};

      template <class... T>
      bool first_non_nil(bool v, T...)
      {
        return v;
      }

      template <class... T>
      bool first_non_nil(nil, T... x)
      {
        return first_non_nil(x...);
      }

    };

    template <class... I>
    struct meta_kernel_context
    {
      typedef std::tuple<I...> image_tuple_t;

      template<class Expr,
	       class Tag = typename proto::tag_of<Expr>::type,
	       class Arg0 = typename Expr::proto_child0>
      struct eval
	: proto::default_eval<Expr, meta_kernel_context, Tag>
      {};

      template <class Expr, int k>
      struct eval<Expr, proto::tag::terminal, ima_expr_tag<k, kernel::point>  >
      {
	typedef typename std::tuple_element<k, image_tuple_t>::type image_t;
	typedef mln_reference(image_t) result_type;

        result_type operator () (Expr&, const meta_kernel_context& ctx) const;
      };

      template <class Expr, int k>
      struct eval<Expr, proto::tag::terminal, ima_expr_tag<k, kernel::neighbor>  >
      {
	typedef typename std::tuple_element<k, image_tuple_t>::type image_t;
	typedef mln_reference(image_t) result_type;

        result_type operator () (Expr&, const meta_kernel_context& ctx) const;
      };

      template <class Expr, class dummy>
      struct eval<Expr, kernel::tag::aggregate, dummy>
      {
	typedef typename proto::result_of::child_c<Expr, 0>::type A0;
	typedef typename proto::result_of::value<A0>::type	  Feature;
	typedef typename proto::result_of::child_c<Expr, 1>::type Arg0;

        typedef typename proto::result_of::eval<Arg0, meta_kernel_context>::type V;
        typedef typename accu::result_of<Feature, typename std::decay<V>::type>::type result_type;

        result_type operator () (Expr&, const meta_kernel_context& ctx) const;
      };

    };

    template <class V, class V2,
              class SubexprTypeList,
              class AccuList>
    struct kernel_context
    {



      kernel_context(AccuList& accus,
                     const V& pval,
                     const V2& nval)
        : m_accus (accus),
          m_pval(pval),
          m_nval(nval)
      {
      }


      template<class Expr,
	       class Tag = typename proto::tag_of<Expr>::type,
	       class Arg0 = typename Expr::proto_child0>
      struct eval
	: proto::default_eval<Expr, const kernel_context, Tag>
      {};


      // Specialization for f(p)
      template <class Expr, int k>
      struct eval<Expr, proto::tag::terminal, ima_expr_tag<k, kernel::point>  >
      {
	typedef typename std::tuple_element<k, V>::type result_type;

	result_type
	operator() (Expr&, const kernel_context& ctx) const
	{
	  return std::get<k>(ctx.m_pval);
	}
      };

      // Specialization for f(n)
      template <class Expr, int k>
      struct eval<Expr, proto::tag::terminal, ima_expr_tag<k, kernel::neighbor>  >
      {
	typedef typename std::tuple_element<k, V>::type result_type;

	result_type
	operator() (Expr&, const kernel_context& ctx) const
	{
	  return std::get<k>(ctx.m_nval);
	}
      };

      // Specialization for terminal A(...) pre-computed
      template <class Expr, int k>
      struct eval<Expr, proto::tag::terminal, kernel::details::aggregate_pl<k>  >
      {
        typedef typename boost::mpl::at_c<SubexprTypeList, k>::type result_type;

	result_type
	operator() (Expr&, const kernel_context& ctx) const
	{
	  return boost::fusion::at_c<k>(ctx.m_accus).to_result();
	}
      };

      // Specialization for A(...) on computation
      // This returns nothing, it only pre-computed an aggragation subexpr
      template <class Expr, class dummy>
      struct eval<Expr, kernel::tag::aggregate, dummy>
      {

	typedef typename proto::result_of::child_c<Expr, 1>::type Arg0;
        typedef typename proto::result_of::eval<Arg0, kernel_context> result_type;

	result_type
	operator() (Expr& subexpr, const kernel_context& ctx) const
	{
	  return proto::eval(proto::child(subexpr), ctx);
	}
      };

    public:
      AccuList& m_accus;
      V         m_pval;
      V2        m_nval;
    };

  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_CONTEXT_HPP
