#ifndef MLN_CORE_IMAGE_CATEGORY_HPP
# define MLN_CORE_IMAGE_CATEGORY_HPP

namespace mln {

  struct forward_image_tag {};
  struct bidirectional_image_tag : forward_image_tag {};
  struct random_access_image_tag : bidirectional_image_tag {};
  struct raw_image_tag : random_access_image_tag {};



} // end of namespace mln

#endif //!MLN_CORE_IMAGE_CATEGORY_HPP
