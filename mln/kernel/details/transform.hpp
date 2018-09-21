#ifndef MLN_KERNEL_DETAILS_TRANSFORM_HPP
# define MLN_KERNEL_DETAILS_TRANSFORM_HPP

# include <boost/proto/proto.hpp>
# include <boost/mpl/apply.hpp>
# include <boost/mpl/int.hpp>
# include <boost/mpl/next.hpp>
# include <boost/proto/functional/fusion/push_back.hpp>
# include <mln/kernel/aggregate.hpp>

namespace mln
{
  namespace kernel
  {
    namespace details
    {
      namespace proto = boost::proto;

      /// \brief The place holder for aggregates
      template <int i>
      struct aggregate_pl : boost::mpl::int_<i>
      {
        friend
        std::ostream&
        operator<<(std::ostream &s, aggregate_pl)
        {
          return s << "Aggr<" << i << ">";
        }
      };


      /// \brief A transform to extract and replace aggregates
      ///
      /// This transform allows to extract any aggregate sub-expression
      /// and to store it in mpl sequence. This aggregate is then replaced
      /// in the AST by aggregate_pl<i> where i is the i nth aggregate.
      ///
      /// \code
      /// auto expr = f(p) + Sum(g(n) - g(p)) + Sum( f(p) )
      /// transform_aggregate trans;
      /// auto texpr = trans(expr)
      template <class start = boost::mpl::int_<0> >
      struct transform_aggregate;

      struct retrieve_aggregate;


      /******************************************/
      /****          Implementation          ****/
      /******************************************/

      struct meta_is_aggregate
      {

        template <class Y, class dummy = void>
        struct apply : boost::mpl::false_ {};

	template <class dummy>
        struct apply< kernel::tag::aggregate, dummy >
          : boost::mpl::true_ {};
      };

      template <class X>
      struct is_aggregate
        : proto::if_< boost::mpl::apply<meta_is_aggregate, proto::tag_of<X> > () >
      {
      };

      template <class start>
      struct meta_tr_aggr : proto::transform< meta_tr_aggr<start> >
      {
        // State is the list of aggragate
        template <class Expr, class State, class Data>
        struct impl : proto::transform_impl<Expr, State, Data>
        {
          typedef typename proto::terminal< aggregate_pl<start::value> >::type result_type;

          result_type
          operator () (typename impl::expr_param,
                       typename impl::state_param,
                       typename impl::data_param)
          {
            return result_type {{}};
          }
        };
      };


      struct count_aggregate_grammar
        : proto::or_<
        proto::when< is_aggregate<proto::_>,
                     boost::mpl::next<proto::_state> ()>,
        proto::when< proto::terminal<proto::_>,
                     proto::_state>,
        proto::otherwise< proto::fold<proto::_, proto::_state, count_aggregate_grammar> >
        >
      {
      };

      template <class Expr>
      struct count_aggregate
      {
        typedef typename std::result_of< count_aggregate_grammar(Expr, boost::mpl::int_<0>) >::type type;
      };

      template <class start, class X>
      struct meta_transform_aggregate
      {
        typedef typename proto::result_of::child<X>::type child0;
        typedef proto::binary_expr<proto::_,
                                   transform_aggregate<start>,
                                   transform_aggregate<
                                     typename std::decay<
                                       typename std::result_of<count_aggregate_grammar(child0, start)>::type
                                       >::type
                                     >
                                   > Transform;

        typedef typename std::result_of<Transform(X)>::type result_type;

        struct type : result_type
        {
          template <typename Expr>
          type(Expr&& x)
            : result_type( Transform() (std::forward<Expr>(x)) )
          {
          }
        };


      };

      template <class start>
      struct transform_aggregate
        : proto::callable, proto::or_<
        proto::when< is_aggregate<proto::_>,
                     meta_tr_aggr<start> >,
        proto::terminal<proto::_>,
        proto::unary_expr<proto::_, transform_aggregate<start> >,
        proto::when<  proto::binary_expr<proto::_, proto::_, proto::_>,
                      meta_transform_aggregate<start, proto::_> (proto::_) >
        >
      {
      };


      struct retrieve_aggregate
        : proto::or_<
        proto::when< is_aggregate<proto::_>,
                     proto::functional::push_back (proto::_state, proto::_)
                     >,
        proto::when< proto::terminal<proto::_>,
                     proto::_state>,
        proto::otherwise< proto::fold<proto::_, proto::_state, retrieve_aggregate> >
        >
      {
      };


    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!MLN_KERNEL_DETAILS_TRANSFORM_HPP
