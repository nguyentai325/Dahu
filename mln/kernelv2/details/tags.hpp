#ifndef MLN_KERNELV2_DETAILS_TAGS_HPP
# define MLN_KERNELV2_DETAILS_TAGS_HPP

# include <iostream>

namespace mln
{

  namespace kernel
  {

    namespace tag
    {

      struct point {};
      struct neighbor {};

      template <int k> struct image_call_n {};
      template <int k> struct image_call_p {};
      struct aggregate {};


      /******************************************/
      /****          Implementation          ****/
      /******************************************/


      template <int k>
      std::ostream&
      operator << (std::ostream& os, image_call_p<k>) {
        return os << "I<" << k << ">(p)";
      }

      template <int k>
      std::ostream&
      operator << (std::ostream& os, image_call_n<k>) {
        return os << "I<" << k << ">(n)";
      }

      template <int k>
      std::ostream&
      operator << (std::ostream& os, aggregate) {
        return os << "A";
      }


    } // end of namespace mln::kernel::tags

  } // end of namespace mln::kernel

} // end of namespace mln

#endif //!MLN_KERNELV2_DETAILS_TAGS_HPP
