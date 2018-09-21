#ifndef MLN_KERNEL_INTRO_HPP
# define MLN_KERNEL_INTRO_HPP

# include <boost/proto/proto.hpp>
# include <boost/mpl/int.hpp>
# include <boost/mpl/max.hpp>
# include <boost/mpl/pair.hpp>
# include <boost/mpl/if.hpp>

# include <mln/core/internal/intseq.hpp>
# include <mln/core/dontcare.hpp>
# include <mln/kernel/types.hpp>

# include <iostream>

/// \file
/// \brief This file provides routines to introspect and transform kernel expression


namespace mln
{

  namespace kernel
  {
    namespace proto = boost::proto;


    /// \brief Provides information about an image usage in a kernel expression.
    ///
    /// Provides information about image usage in a kernel expression. Let \p f
    /// be the \p i nth image in the expression.
    /// It
    /// defines two inner typedefs:
    /// * \p with_p is defined to std::true_type if f(p) appears in the expression
    /// * \p with_n is defined to std::true_type if f(n) appears in the expression
    ///
    /// Exemple:
    /// \code
    /// namespace pl = kernel::placeholders
    /// auto expr = (pl::h(p) = pl::f(pl::p) - kernel::Sum(pl::g(pl::n)));
    /// typedef decltype(expr) E;
    /// std::cout << kernel::get_image_usage<E, 0>::with_p::value; // True
    /// std::cout << kernel::get_image_usage<E, 0>::with_n::value; // False
    /// std::cout << kernel::get_image_usage<E, 1>::with_p::value; // False
    /// std::cout << kernel::get_image_usage<E, 1>::with_n::value; // Tru
    /// std::cout << kernel::get_image_usage<E, 2>::with_p::value; // True
    /// std::cout << kernel::get_image_usage<E, 2>::with_n::value; // False
    /// \endcode
    ///
    /// \tparam
    template <typename Expr, int i>
    struct get_image_usage;


    /// \brief Provides information about a kernel expression.
    ///
    /// Provides information about image usage in a kernel expression. Let \p f
    /// be the \p i nth image in the expression. It defines many typedefs for
    /// each image:
    ///
    /// * \p count is an integral constant for the number of images used in the expression
    /// * \p image<i>::with_p is defined to std::true_type if f(p) appears in the expression
    /// * \p image<i>::with_n is defined to std::true_type if f(n) appears in the expression
    /// Exemple:
    /// \code
    /// namespace pl = kernel::placeholders
    /// auto expr = (pl::h(p) = pl::f(pl::p) - kernel::Sum(pl::g(pl::n)));
    /// typedef decltype(expr) E;
    /// std::cout << kernel::expression_traits<E>::count; // True
    /// std::cout << kernel::get_image_usage<E, 0>::with_n::value; // False
    /// \endcode
    template <typename Expr>
    struct expression_traits;

    /******************************************/
    /****          Implementation          ****/
    /******************************************/


    template <typename Expr, int i>
    struct get_image_usage
    {
    private:
      struct p_grammar
        : proto::or_<
        proto::when<proto::terminal< ima_expr_tag<i, kernel::point> >,
                     boost::mpl::pair<std::true_type, boost::mpl::second<proto::_state> > ()>,
        proto::when<proto::terminal< ima_expr_tag<i, kernel::neighbor> >,
                    boost::mpl::pair<boost::mpl::first<proto::_state>, std::true_type> ()>,
        proto::when<proto::terminal<proto::_>, proto::_state>,
        proto::otherwise<proto::fold<proto::_, proto::_state, p_grammar> >
        >
      {
      };

      typedef typename std::result_of<
        p_grammar(Expr&, boost::mpl::pair<std::false_type, std::false_type>)
        >::type R;

    public:
      typedef typename std::decay<R>::type::first with_p;
      typedef typename std::decay<R>::type::second with_n;
    };


    namespace internal
    {

      struct is_image_expr
      {
        template <typename X>
        struct apply : boost::mpl::false_ {
          //typedef typename X::fuck type;
        };

        template <int i, class P>
        struct apply< ima_expr_tag<i,P> &> : boost::mpl::true_ {
          //typedef typename P::fuck type;
        };

        template <int i, class P>
        struct apply< ima_expr_tag<i,P> > : boost::mpl::true_ {
          //typedef typename P::fuck type;
        };
      };

      template <typename X>
      struct is_image_expr_grammar :
        proto::and_< proto::terminal<X>,
                     proto::if_< boost::mpl::apply<is_image_expr, proto::result_of::value<X> > () >
                     >
      {
      };


      template <class Expr>
      struct get_number_of_images
      {
        struct grammar
          : proto::or_<
          proto::when< is_image_expr_grammar<proto::_>,
                       boost::mpl::max<proto::_state, boost::mpl::next<proto::_value> > ()
                       >,
          proto::when< proto::terminal<proto::_>,
                       proto::_state >,
          proto::otherwise< proto::fold<proto::_, proto::_state,  grammar> >
          >
        {
        };

        typedef typename std::decay<
          typename std::result_of<
            grammar(Expr&, boost::mpl::int_<0>)
            >::type
          >::type type;
      };
    }

    template <typename Expr>
    struct expression_traits
    {
      template <int i>
        using image = get_image_usage<Expr, i>;
      
      typedef typename internal::get_number_of_images<Expr>::type count;
    };




    namespace internal
    {

      template <class Expr, class S>
      struct disp_ker_ima_info;

      template <class Expr, int... I>
      struct disp_ker_ima_info<Expr, intseq<I...> >
      {

        template <int i>
        static void display_()
        {
          std::cout << "Image " << i << ":" << std::endl;
          std::cout << "  use_p:" << expression_traits<Expr>::template image<i>::with_p::value << std::endl;
          std::cout << "  use_n:" << expression_traits<Expr>::template image<i>::with_n::value << std::endl;
        }

        static void display()
        {
          dontcare( (display_<I>(), 0)... );
        }

      };

    }


    template <typename Expr>
    void
    display_kernel(const Expr&)
    {
      typedef typename expression_traits<Expr>::count N;

      std::cout << "The expression is using: " << N::value << " images." << std::endl;

      internal::disp_ker_ima_info<Expr, typename int_list_seq<N::value>::type>::display();

    };


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_INTRO_HPP
