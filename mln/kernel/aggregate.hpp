#ifndef MLN_KERNEL_AGGREGATE_HPP
# define MLN_KERNEL_AGGREGATE_HPP

# include <boost/proto/proto.hpp>

namespace mln
{

  namespace kernel
  {

    namespace tag
    {

      /// \brief The node tag for aggregation nodes.
      struct aggregate {};

    }

    template <class Feature, class Expr>
    using aggregate_expr = typename boost::proto::result_of::make_expr<kernel::tag::aggregate, Feature, Expr>::type;

    template <class Feature>
    struct Aggregate
    {
      template <typename... TParams>
      Aggregate(TParams&&... params)
	: feature(std::forward<TParams>(params) ...)
      {
      }

      template <typename Expr>
      aggregate_expr<Feature, Expr>
      operator() (Expr&& expr) const
      {
        return boost::proto::make_expr<kernel::tag::aggregate, Feature, Expr>(feature, std::forward<Expr>(expr));
      }

      Feature feature;
    };


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_AGGREGATE_HPP
