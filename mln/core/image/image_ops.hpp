#ifndef MLN_CORE_IMAGE_IMAGE_OPS_HPP
# define MLN_CORE_IMAGE_IMAGE_OPS_HPP

# include <mln/core/image/image_expr.hpp>
# include <mln/core/ops.hpp>
# include <boost/utility/enable_if.hpp>
# include <iostream>
namespace mln
{
  /// \defgroup ima_ops Image operators
  /// \{

  /// \brief Check if all elements of image eval to true.
  ///
  /// \param ima The image to check.
  /// \return True iif all elements eval to True.
  ///
  /// \tparam I The type of image. It must be a model of IterableImage.
  /// \pre Image's value type must be convertible to bool.

  template <typename I>
  bool all(const Image<I>& ima);


  /// \brief Check if any element of image eval to true
  ///
  /// \param ima The image to check.
  /// \return True if at least one element evals to True.
  ///
  /// \tparam I The type of image. It must be a model of IterableImage.
  /// \pre Image's value type must be convertible to bool.
  template <typename I>
  bool any(const Image<I>& ima);


  /// \}

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {

    template <typename I, typename Scalar, template<class, class> class Op>
    struct binary_image_scalar_expr_helper
    {
      typedef binary_image_scalar_expr< Op<mln_value(I), Scalar>, I, Scalar> type;
    };

    template <typename I, typename Scalar, template <class, class> class Op>
    struct binary_scalar_image_expr_helper
    {
      typedef binary_scalar_image_expr< Op<Scalar, mln_value(I)>, Scalar, I> type;
    };

  }

# define MLN_GENERATE_LVALUE_UNARY_EXPR(Name, Obj)                      \
  template <typename I>                                                 \
  unary_image_expr< Obj<typename I::value_type>, I>                     \
  Name (I&& ima)                                                        \
  {                                                                     \
    BOOST_CONCEPT_ASSERT(( Image<typename std::remove_reference<I>::type )); \
    return make_unary_image_expr(std::forward<I>(ima), Obj<typename I::reference> ()); \
  }


# define MLN_GENERATE_CONST_UNARY_EXPR(Name, Obj)                       \
  template <typename I>                                                 \
  unary_image_expr< Obj<typename I::value_type>, const I&>                     \
  Name (const Image<I>& ima)                                            \
  {                                                                     \
    return make_unary_image_expr(exact(ima), Obj<typename I::value_type> ()); \
  }

# define MLN_GENERATE_CONST_BINARY_EXPR(Name, Obj)                      \
  template <typename I1, typename I2>                                   \
  binary_image_expr<Obj<mln_value(I1), mln_value(I2)>, I1, I2> \
  Name (const Image<I1>& ima1, const Image<I2>& ima2)                   \
  {                                                                     \
    typedef Obj<typename I1::value_type, typename I2::value_type> O;	\
    return make_binary_image_expr<O, I1, I2>( I1(exact(ima1)), I2(exact(ima2)), O()); \
  }                                                                     \
                                                                        \
                                                                        \
  template <typename I, typename Scalar>                                \
  typename boost::lazy_enable_if_c                                      \
  < !is_a<Scalar, Image>::value,                                        \
    internal::binary_image_scalar_expr_helper<I, Scalar, Obj> >::type	\
  Name (const Image<I>& ima, const Scalar& x)                           \
  {                                                                     \
    typedef Obj<typename I::value_type, Scalar> O;                      \
    return make_binary_image_scalar_expr<O, I, Scalar>(I(exact(ima)), x, O()); \
  }                                                                     \
                                                                        \
  template <typename I, typename Scalar>                                \
  typename boost::lazy_enable_if_c                                      \
  < !is_a<Scalar, Image>::value,                                        \
    internal::binary_scalar_image_expr_helper<I, Scalar, Obj> >::type	\
  Name (const Scalar& x, const Image<I>& ima)                           \
  {                                                                     \
    typedef Obj<Scalar, typename I::value_type> O;                      \
    return make_binary_scalar_image_expr<O, Scalar, I>(x, I(exact(ima)), O()); \
  }



  // Point-wise unary operators
  MLN_GENERATE_CONST_UNARY_EXPR(operator-, negate);
  MLN_GENERATE_CONST_UNARY_EXPR(lnot, logical_not);

  // Point-wise arithmetic operators
  MLN_GENERATE_CONST_BINARY_EXPR(operator+, add);
  MLN_GENERATE_CONST_BINARY_EXPR(operator-, substract);
  MLN_GENERATE_CONST_BINARY_EXPR(operator*, multiplies);
  MLN_GENERATE_CONST_BINARY_EXPR(operator/, devides);
  MLN_GENERATE_CONST_BINARY_EXPR(operator%, modulo);

  // Point-wise relational operators
  MLN_GENERATE_CONST_BINARY_EXPR(operator==, equal_to);
  MLN_GENERATE_CONST_BINARY_EXPR(operator!=, not_equal_to);
  MLN_GENERATE_CONST_BINARY_EXPR(operator<,  less_than);
  MLN_GENERATE_CONST_BINARY_EXPR(operator>,  greater_than);
  MLN_GENERATE_CONST_BINARY_EXPR(operator<=, less_equal);
  MLN_GENERATE_CONST_BINARY_EXPR(operator>=, greater_equal);
  MLN_GENERATE_CONST_BINARY_EXPR(imin, functional::min_t);
  MLN_GENERATE_CONST_BINARY_EXPR(imax, functional::max_t);
  MLN_GENERATE_CONST_BINARY_EXPR(iinf, functional::inf_t);
  MLN_GENERATE_CONST_BINARY_EXPR(isup, functional::sup_t);


  // Point-wise logical operators
  MLN_GENERATE_CONST_BINARY_EXPR(land, logical_and);
  MLN_GENERATE_CONST_BINARY_EXPR(lor, logical_or);



  template <typename I>
  inline
  bool
  all(const Image<I>& ima)
  {
    BOOST_CONCEPT_ASSERT((IterableImage<I>));
    static_assert(std::is_convertible<typename I::value_type, bool>::value,
                  "Image value type must be convertible to bool");

    mln_viter(v, exact(ima));
    mln_forall(v)
      if (!*v) return false;

    return true;
  }

  template <typename I>
  inline
  bool
  any(const Image<I>& ima)
  {
    BOOST_CONCEPT_ASSERT((IterableImage<I>));
    static_assert(std::is_convertible<typename I::value_type, bool>::value,
                  "Image value type must be convertible to bool");

    mln_viter(v, exact(ima));
    mln_forall(v)
      if (*v) return true;

    return false;
  }




} // end of namespace mln



#endif //!MLN_CORE_IMAGE_IMAGE_OPS_HPP
