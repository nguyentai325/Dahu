#ifndef MLN_CORE_IMAGE_IMAGE_HPP
# warning "This should not be included as a standalone."
# include <mln/core/image/image.hpp>
#endif

#ifndef MLN_CORE_CH_VALUE_HPP
# define MLN_CORE_CH_VALUE_HPP

/// \file
/// \brief Meta-Function to transform images.

# define mln_ch_value(I, V) typename image_ch_value<typename std::decay<I>::type, V>::type
# define mln_concrete(I)    typename image_concrete<typename std::decay<I>::type>::type

namespace mln
{

  template <typename Domain>
  struct image_from_domain;


  template <typename V, class Domain, class... TArgs>
  typename image_from_domain<Domain>::template apply<V>::type
  make_image_from_domain(const Domain& dom, TArgs&&... args)
  {
    typedef typename image_from_domain<Domain>::template apply<V>::type OutputImage;
    return OutputImage(dom, std::forward<TArgs>(args)...);
  }


  template <typename Image>
  struct image_concrete :
    image_from_domain<typename Image::domain_type>::template apply< mln_value(Image) >
  {
  };

  // If no specialization get concrete image from domain
  template <typename Image, typename V>
  struct image_ch_value : image_ch_value<typename image_concrete<Image>::type, V>
  {
  };

}

#endif // !MLN_CORE_CH_VALUE_HPP
