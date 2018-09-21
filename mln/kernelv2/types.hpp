#ifndef MLN_KERNELV2_TYPES_HPP
# define MLN_KERNELV2_TYPES_HPP

# include <mln/kernelv2/details/expressions.hpp>

namespace mln
{
  namespace kernel
  {
    namespace proto = boost::proto;

    /**
     * Defines the kernel expression generator available for the user
     *
     */

    /// \brief Current point expression
    struct Point {};

    /// \brief Neighbor point expression
    struct Neighbor {};

    /// \brief image expression
    template <class I, int k> struct Image;

    /// \brief Aggregate expression
    template <class Feature>  struct Aggregate;

    /// \brief Help function to build image expr

    template <int k, class I>
    Image<I, k> make_image_expr(I&& f);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/


    template <class I, int k>
    struct Image
    {
      typedef typename std::remove_reference<I>::type T;

      details::image_call_p_expr<std::reference_wrapper<T>, k>
      operator() (Point)
      {
        return details::image_call_p_expr<std::reference_wrapper<T>, k> { f };
      }

      details::image_call_n_expr<std::reference_wrapper<T>, k>
      operator() (Neighbor)
      {
        return details::image_call_n_expr<std::reference_wrapper<T>, k> { f };
      }

      I f;
    };

    template <class Feature>
    struct Aggregate
    {
      template <typename... TParams>
      Aggregate(TParams&&... params)
	: feature(std::forward<TParams>(params) ...)
      {
      }

      template <typename Expr>
      details::aggregate_expr<Feature, Expr>
      operator() (Expr&& expr) const
      {
        typedef typename proto::result_of::eval<Expr, details::kernel_abstract_context>::type V_;
        typedef typename std::decay<V_>::type V;
        // return details::aggregate_expr<Feature, Expr>::make(accu::make_accumulator(feature, *(V*)NULL),
        //                                                     std::forward<Expr>(expr))

        return proto::make_expr<tag::aggregate, typename accu::accu_of<Feature, V>::type, Expr>
          (accu::make_accumulator(feature, *(V*)NULL),
           std::forward<Expr>(expr));
      }

      Feature feature;
    };


    template <int k, class I>
    Image<I, k>
    make_image_expr(I&& f)
    {
      static_assert(is_a<I, mln::Image>::value,
                    "The first parameter is not an image.");

      return Image<I,k> { std::forward<I>(f) };
    }


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNELV2_TYPES_HPP
