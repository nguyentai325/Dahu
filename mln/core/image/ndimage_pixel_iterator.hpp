#ifndef MLN_CORE_NDIMAGE_PIXEL_ITERATOR_HPP
# define MLN_CORE_NDIMAGE_PIXEL_ITERATOR_HPP


# include <boost/iterator/reverse_iterator.hpp>

# include <mln/core/image/ndimage_pixel.hpp>
# include <mln/core/iterator/fast_reverse_iterator.hpp>
# include <mln/core/image/internal/nested_loop_iterator.hpp>


namespace mln {

  template <typename T, unsigned dim, typename Image>
  using ndimage_pixel_iterator = internal::nested_loop_iterator<
    internal::domain_point_visitor< point<short, dim> >,
    internal::strided_pointer_value_visitor<char, point<short, dim> >,
    ndimage_pixel<T, dim, Image>,
    internal::deref_return_structure_policy>;

template <typename Image, typename Value>
struct ndimage_pixel_range
{
  private:
    enum { ndim = Image::point_type::ndim };
    typedef ndimage_pixel<Value, Image::ndim, Image> pixel_t;

  public:
  typedef internal::nested_loop_iterator<
    internal::origin_point_visitor_forward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim, true>,
    pixel_t,
    internal::deref_return_value_policy> iterator;

  typedef iterator const_iterator;


  typedef internal::nested_loop_iterator<
    internal::origin_point_visitor_backward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim, false>,
    pixel_t,
    internal::deref_return_value_policy> reverse_iterator;

  typedef reverse_iterator const_reverse_iterator;

    ndimage_value_range(Image& ima)
    : ima_(ima)
    {
    }

    iterator iter() const
    {
      Image::point_type pmin = ima_->domain().pmin;
      Image::point_type pmax = ima_->domain().pmax;
      return iterator(pixel_t(ima_),
		      internal::make_point_visitor(pmax - pmin),
		      internal::strided_pointer_value_visitor<ndim>(ima->ptr_, strides.begin()));
    }

    reverse_iterator riter() const
    {
      Image::point_type pmin = ima_->domain().pmin;
      Image::point_type pmax = ima_->domain().pmax;
      return reverse_iterator(pixel_t(ima_),
                              internal::make_point_visitor(pmax - pmin),
                              internal::strided_pointer_value_visitor<ndim>(ima->last_, strides.begin()));
    }

  private:
    Image* ima_;
  };


}


#endif //!MLN_CORE_NDIMAGE_PIXEL_ITERATOR_HPP

