#ifndef MLN_CORE_IMAGE_INTERNAL_REINDEX_HPP
# define MLN_CORE_IMAGE_INTERNAL_REINDEX_HPP


namespace mln
{

  /// \fn void reindex(Image<I>& ima1, const Image<J>& ima2);
  /// \brief Reindex \p ima1 to match indexes of ima2.
  ///
  /// It ensures that after reindexation, \p ima1.domain().pmin and \p
  /// ima2.domain().pmin have the same index. If \p ima1 and/or \p ima2 are not
  /// indexable, the function has no effect.
  ///
  /// Reindexation can be usefull when using view e.g. subimage and still be
  /// able to use indexes.
  ///
  /// \param ima1 The image to be reindexed.
  /// \param ima2 The reference image.
  ///
  ///
  template <typename I, typename J>
  typename std::enable_if<image_traits<I>::indexable::value and image_traits<J>::indexable::value>::type
  reindex(Image<I>& ima1, const Image<J>& ima2);

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename I, typename J>
  inline
  typename std::enable_if<image_traits<I>::indexable::value and image_traits<J>::indexable::value>::type
  reindex(Image<I>& ima, const Image<J>& ima2_)
  {
    const J& ima2 = exact(ima2_);
    auto domain = ima2.domain();
    auto p = domain.iter();
    p.init();
    exact(ima).reindex(ima2.index_of_point(*p));
  }

  template <typename I, typename J>
  inline
  typename std::enable_if<!(image_traits<I>::indexable::value and image_traits<J>::indexable::value)>::type
  reindex(Image<I>&, const Image<J>&)
  {
  }


} // end of namespace mln

#endif //!MLN_CORE_IMAGE_INTERNAL_REINDEX_HPP
