#ifndef MLN_KERNEL_TYPES_HPP
# define MLN_KERNEL_TYPES_HPP

# include <boost/proto/proto.hpp>
# include <boost/mpl/int.hpp>

namespace mln
{
  namespace kernel
  {

    /// \brief Type for point placeholders
    struct point {};

    /// \brief Type for neighbor placeholders
    struct neighbor {};


    template <int I, class P>
    struct ima_expr_tag : boost::mpl::int_<I>
    {
      typedef P arg;

      friend
      std::ostream&
      operator<<(std::ostream &s, ima_expr_tag)
      {
        return s << "Ima<" << I << ">"
                 << "(" << (std::is_same<P,point>::value ? 'p' : 'n') << ")";
      }
    };

    /// \brief Node type for a pixel value access through a point or a neighbor
    template <int I, class P>
    using ima_expr = typename boost::proto::terminal< ima_expr_tag<I, P> >::type;


  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNEL_TYPES_HPP
