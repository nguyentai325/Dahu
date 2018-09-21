#ifndef MLN_CORE_IMAGE_IMAGE_HPP
# warning "You should not include this file directly but <mln/core/image/image.hpp> instead"
# include <mln/core/image/image.hpp>
#endif

#ifndef MLN_CORE_IMAGE_IMAGE_EXPR_HPP
# define MLN_CORE_IMAGE_IMAGE_EXPR_HPP

# include <type_traits>

# include <mln/core/image/morphers/zip_image.hpp>
# include <mln/core/image/morphers/transformed_image.hpp>

namespace mln
{

  namespace internal
  {
    template <typename BinaryFunction, class U, class V>
    struct func_call_from_tupleargs;
  }

  template <typename UnaryFunction, typename Image>
  using unary_image_expr = internal::transformed_image<Image, UnaryFunction, false>;

  template <typename BinaryFunction, typename Image1, typename Image2>
  using binary_image_expr = internal::transformed_image<
    zip_image<Image1, Image2>,
    internal::func_call_from_tupleargs<BinaryFunction, mln_reference(Image1), mln_reference(Image2)> >;

  template <typename BinaryFunction, typename Image, typename Scalar>
  using binary_image_scalar_expr = internal::transformed_image<
    Image,
    decltype( std::bind(std::declval<BinaryFunction>(), std::placeholders::_1, std::declval<Scalar>()) )
    >;

  template <typename BinaryFunction, typename Scalar, typename Image>
  using binary_scalar_image_expr = internal::transformed_image<
    Image,
    decltype( std::bind(std::declval<BinaryFunction>(), std::declval<Scalar>(), std::placeholders::_1) )
    >;

  // template <typename UnaryFunction, typename I>
  // unary_image_expr<UnaryFunction, I>
  // make_unary_image_expr(I&& ima, UnaryFunction f);




  /*
  template <typename I>
  I& eval(Image<I>& ima);

  template <typename I>
  const I& eval(const Image<I>& ima);

  template <typename Expr>
  typename eval_type<Expr>::type
  eval(const image_expr<Expr>& ima)
  */


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {

    template <typename Function, class U, class V>
    struct func_call_from_tupleargs
    {
      func_call_from_tupleargs() = default;

      func_call_from_tupleargs(const Function& f_)
      : f(f_)
      {
      }

      typename std::result_of<Function(U, V)>::type
      operator() (std::tuple<U, V> x) const
      {
        return f(std::get<0>(x), std::get<1>(x));
      }

    private:
      Function f;
    };

  }


  /******************************************/
  /****         Helper functions         ****/
  /******************************************/

  template <typename UnaryFunction, typename I>
  unary_image_expr<UnaryFunction, I>
  make_unary_image_expr(I&& ima, const UnaryFunction& f)
  {
    BOOST_CONCEPT_ASSERT((Image<typename std::decay<I>::type>));
    return unary_image_expr<UnaryFunction, I>(std::forward<I>(ima), f);
  }

  template <typename BinaryFunction, typename I1, typename I2>
  binary_image_expr<BinaryFunction, I1, I2>
  make_binary_image_expr(I1&& ima1, I2&& ima2, const BinaryFunction& f)
  {
    BOOST_CONCEPT_ASSERT((Image<typename std::decay<I1>::type>));
    BOOST_CONCEPT_ASSERT((Image<typename std::decay<I2>::type>));
    auto z = imzip(std::forward<I1>(ima1), std::forward<I2>(ima2));
    internal::func_call_from_tupleargs<BinaryFunction, mln_reference(I1), mln_reference(I2)> fun(f);
    return binary_image_expr<BinaryFunction, I1, I2>(std::move(z), fun);
  }


  template <typename BinaryFunction, typename I, typename Scalar>
  binary_image_scalar_expr<BinaryFunction, I, Scalar>
  make_binary_image_scalar_expr(I&& ima, const Scalar& x, const BinaryFunction& f)
  {
    BOOST_CONCEPT_ASSERT((Image<typename std::decay<I>::type>));
    auto fun = std::bind(f, std::placeholders::_1, x);
    return binary_image_scalar_expr<BinaryFunction, I, Scalar>(std::forward<I>(ima), fun);
  }

  template <typename BinaryFunction, typename Scalar, typename I>
  binary_scalar_image_expr<BinaryFunction, Scalar, I>
  make_binary_scalar_image_expr(const Scalar& x, I&& ima, const BinaryFunction& f)
  {
    BOOST_CONCEPT_ASSERT((Image<typename std::decay<I>::type>));
    auto fun = std::bind(f, x, std::placeholders::_1);
    return binary_scalar_image_expr<BinaryFunction, Scalar, I>(std::forward<I>(ima), fun);
  }

}

#endif // !MLN_CORE_IMAGE_IMAGE_EXPR_HPP
