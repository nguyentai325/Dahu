#ifndef MLN_CORE_IMAGE_MORPHERS_CASTED_IMAGE_HPP
# define MLN_CORE_IMAGE_MORPHERS_CASTED_IMAGE_HPP

# include <mln/core/image/morphers/transformed_image.hpp>

namespace mln
{
  // Exposition only
  namespace internal
  {
    template <typename Vout>
    struct cast_to;
  }
  // End

  template <class I, class V>
  using casted_image = internal::transformed_image<I, internal::cast_to<V>, false>;

  template <class V, class I>
  casted_image<const I&, V>
  imcast(const Image<I>& ima);

  template <class V, class I>
  casted_image<I&, V>
  imcast(Image<I>& ima);

  template <class V, class I>
  casted_image<I&&, V>
  imcast(Image<I>&& ima);

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {
    template <typename V>
    struct cast_to
    {
      template <typename T>
      V
      operator() (const T& v) const
      {
	return static_cast<V>(v);
      }
    };
  }

  template <class V, class I>
  casted_image<const I&, V>
  imcast(const Image<I>& ima)
  {
    return casted_image<const I&, V> (exact(ima), internal::cast_to<V> ());
  }

  template <class V, class I>
  casted_image<I&, V>
  imcast(Image<I>& ima)
  {
    return casted_image<I&, V> (exact(ima), internal::cast_to<V> ());
  }

  template <class V, class I>
  casted_image<I&&, V>
  imcast(Image<I>&& ima)
  {
    return casted_image<I&&, V> (move_exact(ima), internal::cast_to<V> ());
  }

}

#endif // ! MLN_CORE_IMAGE_MORPHERS_CASTED_IMAGE_HPP
