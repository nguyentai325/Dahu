#ifndef MLN_CORE_SLIDING_IMAGE_ITERATOR_HPP
# define MLN_CORE_SLIDING_IMAGE_ITERATOR_HPP

# include <mln/core/pixel.hpp>
# include <boost/iterator/iterator_facade.hpp>
# include <mln/core/std/array.hpp>
# include <mln/core/image_traits.hpp>

namespace mln {


  template <typename Image, typename SiteSet, typename image_tag = typename image_traits<Image>::category>
  struct sliding_image_value_iterator;


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <typename Image, size_t N>
  struct sliding_image_value_iterator<Image, mln::array<typename Image::point_type, N>, raw_image_tag>
    : iterator_base< sliding_image_value_iterator<Image, mln::array<typename Image::point_type, N>, raw_image_tag>,
		     typename Image::value_type,
		     typename image_reference<Image>::type>
  {
  private:
    typedef std::array<typename Image::difference_type, N> offset_t;
    typedef typename image_pointer<Image>::type		ptr_t;

  public:
    typedef typename image_reference<Image>::type	reference;

    sliding_image_value_iterator()
    {
    }

    sliding_image_value_iterator(ptr_t ptr, const offset_t& offsets)
      : offsets_ (offsets), ptr_ ((char*)ptr)
    {
    }


    void init() { i_ = 0; }
    void next() { ++i_; }
    bool finished() const { return i_ == N; }
    reference dereference() const { return *(ptr_t)(ptr_ + offsets_[i_]); }

  private:
    offset_t offsets_;
    char*    ptr_;
    int      i_;
  };

}

#endif //!MLN_CORE_SLIDING_IMAGE_ITERATOR_HPP
