#ifndef MLN_MORPHO_CONCEPT_SE_HPP
# define MLN_MORPHO_CONCEPT_SE_HPP

# include <mln/core/concept/object.hpp>

namespace mln
{
  namespace morpho
  {

  template <class SE>
  using StructuringElement = Neighborhood<SE>;


  } // end of namespace mln::morpho

} // end of namespace mln

#endif //!MLN_MORPHO_CONCEPT_SE_HPP
