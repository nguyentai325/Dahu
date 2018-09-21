#ifndef MLN_KERNELV2_DEBUG_HPP
# define MLN_KERNELV2_DEBUG_HPP

# include <mln/kernelv2/kernel.hpp>
# include <mln/kernelv2/details/expressions_traits.hpp>
# include <mln/core/dontcare.hpp>

namespace mln
{
  namespace kernel
  {

    template <class Expr>
    void
    debug(Expr&& expr);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/
    namespace internal
    {

      template <int k, class Expr>
      void debug_image(Expr&& x)
      {
        typedef details::image_usage_traits<Expr,k> traits;
        typedef typename traits::type I;

        std::cout << "Image " << k << " {\n"
                  << "  adress: " << &(details::get_image<k>(std::forward<Expr>(x))) << "\n"
                  << "  type: " << typeid(typename traits::type).name()
                  << (std::is_reference<I>::value ? "(ref)" : "") << "\n"
                  << "  point: " << traits::point::value << "\n"
                  << "  neighbor: " << traits::neighbor::value << "\n"
                  << "}" << std::endl;
      }

      template <class Expr, int... k >
      void debug_all_image(Expr&& expr, intseq<k...>)
      {
        dontcare((debug_image<k>(std::forward<Expr>(expr)), (NULL))...);
      }

    }

    template <class Expr>
    void
    debug(Expr&& expr)
    {
      proto::display_expr(std::forward<Expr>(expr));

      typedef typename details::expression_traits<Expr>::number_of_images n;

      std::cout << "Number of images: "
                << n::value
                << std::endl;

      internal::debug_all_image(expr, typename int_list_seq<n::value>::type ());
    }

  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNELV2_DEBUG_HPP
