#ifndef MLN_CORE_LITERAL_VECTORIAL_HPP
# define MLN_CORE_LITERAL_VECTORIAL_HPP

namespace mln
{
  namespace literal
  {

  struct zero_t
  {
    constexpr operator int () const
    {
      return 0;
    }
  };

  struct one_t
  {
    constexpr operator int () const
    {
      return 1;
    }
  };


    static const zero_t zero = {};
    static const one_t  one = {};

  } // end of namespace mln::literal

} // end of namespace mln

#endif //!MLN_CORE_LITERAL_VECTORIAL_HPP
