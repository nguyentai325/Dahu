#ifndef MLN_CORE_IMAGE_SUB_IMAGE_SPE_HPP
# define MLN_CORE_IMAGE_SUB_IMAGE_SPE_HPP

# include <type_traits>

namespace mln
{

  // Specialization for ndimage | box
  // fwd decl
  template <typename, unsigned, typename> struct ndimage_base;

  template <typename T, unsigned dim, typename E, typename Domain>
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, const E>::type
  make_subimage(const ndimage_base<T, dim, E>&,
		const Domain& domain);

  template <typename T, unsigned dim, typename E, typename Domain>
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, E>::type
  make_subimage(ndimage_base<T, dim, E>&,
		const Domain& domain);

  template <typename T, unsigned dim, typename E, typename Domain>
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, E>::type
  make_subimage(ndimage_base<T, dim, E>&&,
		const Domain& domain);


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename T, unsigned dim, typename E, typename Domain>
  inline
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, E>::type
  make_subimage(ndimage_base<T, dim, E>& image,
                const Domain& domain)
  {
    E other(exact(image));
    other.domain_ = domain;
    other.border_ = image.border_; // FIXME
    other.ptr_ = (char*) &image(domain.pmin);
    other.last_ = (char*) &image(domain.pmax - 1);
    other.m_ptr_origin = image.m_ptr_origin;
    other.m_index_first = image.index_of_point(domain.pmin);
    other.m_index_last = image.index_of_point(domain.pmax);
    return other;
  }

  template <typename T, unsigned dim, typename E, typename Domain>
  inline
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, const E>::type
  make_subimage(const ndimage_base<T, dim, E>& image,
                const Domain& domain)
  {
    return make_subimage(const_cast< ndimage_base<T, dim, E>& >(image), domain);
  }


  template <typename T, unsigned dim, typename E, typename Domain>
  inline
  typename std::enable_if< std::is_convertible<Domain, typename ndimage_base<T, dim, E>::domain_type>::value, E>::type
  make_subimage(ndimage_base<T, dim, E>&& image,
                const Domain& domain)
  {
    return make_subimage(image, domain);
  }


} // end of namespace mln

#endif //!MLN_CORE_IMAGE_SUB_IMAGE_SPE_HPP
