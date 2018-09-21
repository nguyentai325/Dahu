#ifndef MLN_CORE_SLIDING_OBJECT_HXX
# define MLN_CORE_SLIDING_OBJECT_HXX

namespace mln {
  namespace internal {

    template <typename Image, int N>
    struct sliding_object_base<Image, (typename Image::point_type)[N], raw_access_image_tag>
    {
      // As an image
      typedef typename Image::point_type        point_type;
      typedef typename Image::value_type        value_type;

      typedef typename Image::difference_type   difference_type;

      // Constructor
      sliding_object_base(Image& ima, point_type (&dpoints)[N]);

      // iterators
      typedef sliding_object_pixel_iterator     



      // For this image only:
      void center(point_type p) { p_ = p; }

    private:
      Image&                ima_;
      point_type&           dpoints[N];
      difference_type       offsets_[N];

      point_type            p_;
    };


  } // end of namespace mln::internal
} // end of namespace mln




#endif //!MLN_CORE_SLIDING_OBJECT_HXX
