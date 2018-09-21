#ifndef MLN_NDIMAGE_RANGE_HH
# define MLN_NDIMAGE_RANGE_HH

# include <mln/core/image/internal/nested_loop_iterator.hpp>
# include <mln/core/image/ndimage_pixel.hpp>
# include <mln/core/stride_utils.hpp>

namespace mln
{

  template <typename Image, typename Value>
  struct ndimage_value_range
  {
  private:
    enum { ndim = Image::point_type::ndim };
    typedef typename exact_type<Image>::type E;
    typedef ndimage_pixel<Value, Image::ndim, E> pixel_t;

    template <typename J, typename V>
    friend struct ndimage_value_range;

  public:

    typedef internal::nested_loop_iterator<
    internal::origin_point_visitor_forward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim>,
    internal::no_op_visitor,
    pixel_t, internal::deref_return_value_policy> iterator;

  typedef internal::nested_loop_iterator<
    internal::origin_point_visitor_backward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim>,
    internal::no_op_visitor,
    pixel_t, internal::deref_return_value_policy> reverse_iterator;

    typedef iterator const_iterator;
    typedef reverse_iterator const_reverse_iterator;


    ndimage_value_range(Image& ima)
    : ima_(&ima)
    {
    }

    // interop
    ndimage_value_range(const typename Image::value_range& other)
      : ima_(other.ima_)
    {
    }


    const_iterator iter() const
    {
      typename Image::point_type pmin = ima_->domain().pmin;
      typename Image::point_type pmax = ima_->domain().pmax;
      typename Image::point_type shp = pmax - pmin;

      std::array<std::ptrdiff_t, ndim> delta_byte_strides;
      delta_byte_strides[ndim-1]  = ima_->strides_[ndim-1];

      for (unsigned i = 0; i < ndim-1; ++i)
        delta_byte_strides[i] = ima_->strides_[i] - ima_->strides_[i+1] * (shp[i+1] - 1);

      return const_iterator(pixel_t(exact(ima_)),
			    internal::make_point_visitor_forward(shp),
			    internal::strided_pointer_value_visitor<ndim>(ima_->ptr_, delta_byte_strides),
                            internal::no_op_visitor ()
			    );
    }

    const_reverse_iterator riter() const
    {
      typename Image::point_type pmin = ima_->domain().pmin;
      typename Image::point_type pmax = ima_->domain().pmax;
      typename Image::point_type shp = pmax - pmin;

      std::array<std::ptrdiff_t, ndim> delta_byte_strides;
      delta_byte_strides[ndim-1]  = ima_->strides[ndim-1];

      for (unsigned i = 0; i < ndim-1; ++i)
        delta_byte_strides[i] = -(ima_->strides_[i] - ima_->strides[i+1] * (shp[i+1] - 1));

      return const_reverse_iterator(pixel_t(exact(ima_)),
				    internal::make_point_visitor_forward(shp),
				    internal::strided_pointer_value_visitor<ndim>(ima_->last_, delta_byte_strides),
                                    internal::no_op_visitor ()
				    );
    }

  private:
    Image* ima_;
  };


  template <typename Image, typename Value>
  struct ndimage_pixel_range
  {
  private:
    enum { ndim = Image::point_type::ndim };
    typedef typename exact_type<Image>::type E;
    typedef ndimage_pixel<Value, Image::ndim, E> pixel_t;

    template <typename I, typename V>
    friend struct ndimage_pixel_range;

  public:
    typedef internal::nested_loop_iterator<
    internal::domain_point_visitor_forward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim>,
    internal::strided_index_visitor<ndim>,
    pixel_t, internal::deref_return_structure_policy> iterator;

    typedef iterator const_iterator;

  typedef internal::nested_loop_iterator<
    internal::domain_point_visitor_backward< typename Image::point_type >,
    internal::strided_pointer_value_visitor<ndim>,
    internal::strided_index_visitor<ndim>,
    pixel_t, internal::deref_return_structure_policy> reverse_iterator;

    typedef reverse_iterator const_reverse_iterator;


    ndimage_pixel_range(Image& ima)
    : ima_(&ima)
    {
    }

    ndimage_pixel_range(const typename Image::pixel_range& other)
      : ima_(other.ima_)
    {
    }

    iterator iter() const
    {
      typename Image::point_type pmin = ima_->domain().pmin;
      typename Image::point_type pmax = ima_->domain().pmax;
      typename Image::point_type shp = pmax - pmin;
      //typename Image::point_type bshp = shp + 2*ima_->border_;

      std::array<std::ptrdiff_t, ndim> delta_byte_strides;
      std::array<std::ptrdiff_t, ndim> delta_index_strides;

      compute_delta_strides<ndim>(& ima_->m_index_strides[0], &shp[0], &delta_index_strides[0]);
      compute_delta_strides<ndim>(& ima_->strides_[0],        &shp[0], &delta_byte_strides[0]);
      // std::cout << ima_->m_index_strides[0] << std::endl;
      // std::cout << ima_->m_index_strides[1] << std::endl;
      // std::cout << delta_index_strides[0] << std::endl;
      // std::cout << delta_index_strides[1] << std::endl;

      return iterator(pixel_t(exact(ima_)),
		      internal::make_point_visitor_forward(pmin, pmax),
		      internal::strided_pointer_value_visitor<ndim>(ima_->ptr_, delta_byte_strides),
		      internal::strided_index_visitor<ndim>(ima_->m_index_first, delta_index_strides)
		      );
    }

    reverse_iterator riter() const
    {
      typename Image::point_type pmin = ima_->domain().pmin;
      typename Image::point_type pmax = ima_->domain().pmax;
      typename Image::point_type shp = pmax - pmin;
      //typename Image::point_type bshp = shp + 2*ima_->border_;

      std::array<std::ptrdiff_t, ndim> delta_byte_strides;
      std::array<std::ptrdiff_t, ndim> delta_index_strides;
      compute_negative_delta_strides<ndim>(&ima_->m_index_strides[0], &shp[0], &delta_index_strides[0]);
      compute_negative_delta_strides<ndim>(&ima_->strides_[0],        &shp[0], &delta_byte_strides[0]);

      return reverse_iterator(pixel_t(exact(ima_)),
			      internal::make_point_visitor_backward(pmin, pmax),
			      internal::strided_pointer_value_visitor<ndim>(ima_->last_, delta_byte_strides),
			      internal::strided_index_visitor<ndim>(ima_->m_index_last, delta_index_strides)
			      );
    }
  private:
    Image* ima_;
  };
}
#endif
