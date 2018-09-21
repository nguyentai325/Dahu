#ifndef MLN_KERNELV2_EXECUTE_HPP
# define MLN_KERNELV2_EXECUTE_HPP

# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/kernelv2/details/expressions_traits.hpp>
# include <mln/kernelv2/details/inner_context.hpp>
# include <mln/kernelv2/details/eval_context.hpp>
# include <mln/kernelv2/details/pixter_wrapper.hpp>

namespace mln
{
  namespace kernel
  {
    ///
    ///
    /// Note that the expression is copied thus if the expression
    /// involves very expensive copy types, use reference wrappers.
    /// You shouls also do so with side-effect types.
    template <class Expr, class N>
    void execute(Expr x, const Neighborhood<N>& nbh);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    template <class Expr, class N>
    void execute(Expr x, const Neighborhood<N>& nbh_)
    {
      const N& nbh = exact(nbh_);

      typedef typename details::image_list_traits<Expr&>::zip_image_type Z_1;
      typedef typename details::image_used_by_neighbor_list_traits<Expr&>::zip_image_type Z_2;


      Z_1 f = details::get_image_list2(x);
      Z_2 g = details::get_image_used_by_neighbor_list2(x);

      typedef mln_reference(Z_1) V1;
      typedef mln_reference(Z_2) V2;

      // Outer iterator is rebinded to a fake image
      // so that the neighborhood won't iterate on an image if the
      // final value is not used.
      mln_pixter(px, f);
      details::wrap_pixter<Z_2, decltype(px)> px__(g, &px);

      // Inner iterator
      mln_iter(nx, nbh(px__));

      mln_forall(px)
      {
        { // Call init on every aggregate
          details::inner_context_init ctx;
          proto::eval(x, ctx);
        }

        mln_forall(nx)
        {
          details::inner_context_take<V1,V2> ctx = {px->val(), nx->val()};
          proto::eval(x, ctx);
        }

        details::eval_context<V1,dontcare_t> ctx = {px->val(), dontcare};
        proto::eval(x, ctx);
      }
    }


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNELV2_EXECUTE_HPP
