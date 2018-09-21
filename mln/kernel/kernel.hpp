#ifndef MLN_KERNEL_HPP
# define MLN_KERNEL_HPP

# include <boost/proto/proto.hpp>
# include <mln/core/image/image.hpp>
# include <mln/kernel/execute.hpp>

namespace mln
{

  namespace kernel
  {

    namespace proto = boost::proto;

    // using Max_t = typename proto::terminal< mln::accu::features::max<> >::type;
    // static const Max_t Max = {{}};


    template <int I>
    struct image
    {
      ima_expr<I, kernel::point>
      operator() (kernel::point) const
      {
	return ima_expr<I, kernel::point> {{}};
      }


      ima_expr<I, kernel::neighbor>
      operator() (kernel::neighbor) const
      {
	return ima_expr<I, kernel::neighbor> {{}};
      };
    };

    /*********************************/
    /****    Placeholders        *****/
    /*********************************/

    namespace placeholders
    {
      static const kernel::image<0> f {};
      static const kernel::image<1> g {};
      static const kernel::image<2> h {};
      static const kernel::image<3> out {};

      static const kernel::point p {};
      static const kernel::neighbor n {};
    };


    /************************************/
    /** Default kernel Implementation  **/
    /************************************/

  }

}

#endif // ! MLN_KERNEL_HPP
