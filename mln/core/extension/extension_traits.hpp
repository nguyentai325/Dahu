#ifndef MLN_CORE_EXTENSION_EXTENSION_TRAITS_HPP
# define MLN_CORE_EXTENSION_EXTENSION_TRAITS_HPP

namespace mln
{
  namespace extension
  {

/******************************************/
/****               Tags               ****/
/******************************************/

    struct extension_tag {};
    struct none_extension_tag : extension_tag {};
    struct custom_extension_tag : extension_tag {};
    struct border_extension_tag : custom_extension_tag {};
    struct value_extension_tag : custom_extension_tag {};
    struct image_extension_tag : custom_extension_tag {};


  } // end of namespace mln::extension

/******************************************/
/****              Traits              ****/
/******************************************/

    template <typename E>
    struct extension_traits
    {
      typedef typename E::support_fill           support_fill;
      typedef typename E::support_mirror         support_mirror;
      typedef typename E::support_periodize      support_periodize;
    };


} // end of namespace mln

#endif //!MLN_CORE_EXTENSION_EXTENSION_TRAITS_HPP
