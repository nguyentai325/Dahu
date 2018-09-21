#ifndef MLN_CORE_ALGORITHM_EQUAL_HPP
# define MLN_CORE_ALGORITHM_EQUAL_HPP

# include <mln/core/image/image.hpp>

namespace mln
{

  /// \brief Compare if two image are equals.
  ///
  /// Two image are said equal if their domains
  /// are equals and have the same values.
  ///
  /// \param ima1 First image
  /// \param ima2 Second image
  ///
  /// \tparam I must model a Readable Forward Image
  /// \tparam J must model a Readable Forward Image
  ///
  /// \return True if image are equals.
  ///
  template <typename I, typename J>
  bool
  equal(const Image<I>& ima1, const Image<J>& ima2);


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename I, typename J>
  inline
  bool
  equal(const Image<I>& ima1, const Image<J>& ima2)
  {
    mln_pixter(px1, exact(ima1));
    mln_pixter(px2, exact(ima2));

    mln_forall(px1, px2)
      if (px1->point() != px2->point() or
          px1->val() != px2->val())
        return false;

    return true;
  };


} // end of namespace mln

#endif //!MLN_CORE_ALGORITHM_EQUAL_HPP
