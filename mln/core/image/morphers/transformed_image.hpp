// Transformed image is core morpher, include image.hpp to ensure dependencies
#include <mln/core/image/image.hpp>

#ifndef MLN_CORE_IMAGE_MORPHERS_TRANSFORMED_IMAGE_HPP
# define MLN_CORE_IMAGE_MORPHERS_TRANSFORMED_IMAGE_HPP

// FIXME: pixel do not preserve properties

# include <mln/core/range/transform.hpp>
# include <mln/core/pixel_utility.hpp>

namespace mln
{
  // Exposition only
  namespace internal
  {
    template <typename I, class UnaryFunction, bool use_pix = use_pix_helper<I, UnaryFunction>::type::value >
    struct transformed_image;
  }
  // End

  template <typename I, class UnaryFunction>
  using transformed_image = internal::transformed_image<I, UnaryFunction>;

  template <typename I, class UnaryFunction>
  transformed_image<I, UnaryFunction>
  imtransform(Image<I>&& ima, const UnaryFunction& f);

  template <typename I, class UnaryFunction>
  transformed_image<const I&, UnaryFunction>
  imtransform(const Image<I>& ima, const UnaryFunction& f);

  template <typename I, class UnaryFunction>
  transformed_image<I&, UnaryFunction>
  imtransform(Image<I>& ima, const UnaryFunction& f);
}

  /******************************************/
  /****           HELPER MACROS          ****/
  /******************************************/

  // template <typename I>                                                 
  // transformed_image<I, BOOST_PP_TUPLE_REM () F, false>                  
  // NAME##_helper(I&& ima, BOOST_PP_TUPLE_REM () TYPE *)                  
  // {                                                                     
  //   BOOST_PP_TUPLE_REM () F f;                                          
  //   return transformed_image<I, BOOST_PP_TUPLE_REM () F, false>(std::forward<I>(ima), f); 
  // }                                                                     



# define  MLN_INTERNAL_IMAGE_LVALUE_OPERATOR_TEMPLATE_2(NAME, TYPE, F)  \
 namespace mln { namespace internal {                                  \
                                                                        \
      template <typename T>                                             \
      struct NAME##_helper_t<BOOST_PP_TUPLE_REM () TYPE>                \
      {                                                                 \
        typedef BOOST_PP_TUPLE_REM () TYPE          type;               \
        typedef BOOST_PP_TUPLE_REM () F             f_type;             \
                                                                        \
        template <class I>                                              \
        struct apply                                                    \
        {                                                               \
          typedef transformed_image<I, f_type, false> result_type;      \
          result_type                                                   \
            operator() (I f) const                                      \
          {                                                             \
            return result_type(std::forward<I>(f), f_type());           \
          }                                                             \
        };                                                              \
      };                                                                \
   }}

# define  MLN_INTERNAL_IMAGE_LVALUE_OPERATOR_2(NAME, TYPE, F)           \
  namespace mln { namespace internal {                                  \
      template<>                                                        \
      struct NAME##_helper_t<BOOST_PP_TUPLE_REM () TYPE>                \
      {                                                                 \
        typedef BOOST_PP_TUPLE_REM () TYPE          type;               \
        typedef BOOST_PP_TUPLE_REM () F             f_type;             \
                                                                        \
        template <class I>                                              \
        struct apply                                                    \
        {                                                               \
          typedef transformed_image<I, f_type, false> result_type;      \
          result_type                                                   \
            operator() (I f) const                                      \
          {                                                             \
            return result_type(std::forward<I>(f), f_type());           \
          }                                                             \
        };                                                              \
      };                                                                \
    }}


# define MLN_INTERNAL_IMAGE_DECLARE_OPERATOR(NAME)                      \
                                                                        \
  namespace mln {                                                       \
                                                                        \
    namespace internal {                                                \
      template <class T>                                                \
      struct NAME##_helper_t;                                           \
    }                                                                   \
                                                                        \
    template <typename I>                                               \
    typename internal::NAME##_helper_t<mln_value(I)>::template apply<I&&>::result_type \
    NAME(Image<I>&& ima)                                                \
    {									\
      return typename internal::NAME##_helper_t<mln_value(I)>::template apply<I&&> () (move_exact<I>(ima)); \
    }									\
                                                                        \
    template <typename I>                                               \
    typename internal::NAME##_helper_t<mln_value(I)>::template apply<I&>::result_type \
    NAME(Image<I>& ima)                                                 \
    {									\
      return typename internal::NAME##_helper_t<mln_value(I)>::template apply<I&> () (exact(ima)); \
    }									\
                                                                        \
    template <typename I>                                               \
    typename internal::NAME##_helper_t<mln_value(I)>::template apply<const I&>::result_type \
    NAME(const Image<I>& ima)                                           \
    {									\
      return typename internal::NAME##_helper_t<mln_value(I)>::template apply<const I&> () (exact(ima)); \
    }                                                                   \
  }



# define  MLN_DECLARE_IMAGE_LVALUE_OPERATOR_OVERLOAD(NAME, TYPE, F)     \
  MLN_INTERNAL_IMAGE_LVALUE_OPERATOR_2(NAME, TYPE, F)


# define  MLN_DECLARE_IMAGE_LVALUE_OPERATOR(NAME, TYPE, F)              \
  MLN_INTERNAL_IMAGE_DECLARE_OPERATOR(NAME);                            \
  MLN_INTERNAL_IMAGE_LVALUE_OPERATOR_2(NAME, TYPE, F)



# define  MLN_DECLARE_IMAGE_LVALUE_OPERATOR_TEMPLATE(NAME, TYPE, F)     \
  MLN_INTERNAL_IMAGE_DECLARE_OPERATOR(NAME)                             \
  MLN_INTERNAL_IMAGE_LVALUE_OPERATOR_TEMPLATE_2(NAME, TYPE, F)



  /******************************************/
  /****              Traits              ****/
  /******************************************/
namespace mln
{

  template <typename I, class UnaryFunction, bool b>
  struct image_traits< internal::transformed_image<I, UnaryFunction, b> >
  {
  private:
    typedef typename std::remove_reference<I>::type image_t;

  public:
    typedef typename std::common_type<typename image_traits<image_t>::category,
                                      random_access_image_tag>::type category;
    typedef std::false_type				 concrete;
    typedef typename image_traits<image_t>::accessible   accessible;
    typedef typename image_traits<image_t>::indexable    indexable;
    typedef mln::extension::none_extension_tag		 extension; // FIXME
  };

  template <class I, class UnaryFunction, bool b>
  struct image_concrete< internal::transformed_image<I, UnaryFunction, b> >
  {
    typedef typename std::result_of<UnaryFunction(mln_reference(I))>::type __vtype;
    typedef typename std::decay<__vtype>::type                             value_type;

    typedef mln_ch_value(typename std::remove_reference<I>::type, value_type) type;
  };

  namespace internal
  {
    template <class I, class UnaryFunction, bool b>
    struct image_init_from< transformed_image<I, UnaryFunction, b> >
    {
      typedef typename image_init_from<typename std::decay<I>::type>::type type;
    };
  }

  /******************************************/
  /****          Implementation          ****/
  /******************************************/


  namespace internal
  {
    // Helper for clang, it refuses to hold mutable reference
    template <typename I, class Enable = typename std::is_reference<I>::type>
    struct transformed_image_helper;

    template <typename I>
    struct transformed_image_helper<I, std::true_type>
    {
      I m_ima;
    };

    template <typename I>
    struct transformed_image_helper<I, std::false_type>
    {
      mutable I m_ima;
    };

    // With a value transformation: V -> W
    template <typename I, class UnaryFunction>
    struct transformed_image<I, UnaryFunction, false>
      : morpher_base< transformed_image<I, UnaryFunction, false>,
                      I, mln_point(I),
                      typename std::remove_reference
                      <typename std::result_of<UnaryFunction(mln_reference(I))>::type>::type >,
        transformed_image_helper<I>
    {
    private:
      typedef typename std::remove_reference<I>::type image_t;
      typedef transformed_image<I, UnaryFunction, false> this_t;
      friend struct mln::morpher_core_access;
      typedef typename std::result_of<UnaryFunction(mln_reference(I))>::type    _result_type;

    public:
      // If F(x) -> T&& then the return type is T (because it would yield a dangling rvalue ref to x)
      // reference is thus T&, const T& or T
      typedef typename std::conditional<std::is_rvalue_reference<_result_type>::value,
                                        typename std::remove_reference<_result_type>::type,
                                        _result_type>::type                      reference;

      typedef typename std::decay<reference>::type				value_type;


      // if reference is T& or T const &, the const_reference is T const &
      // if reference is T, the const_reference is still T
      typedef typename std::conditional<
        std::is_reference<reference>::value,
        typename std::add_lvalue_reference<typename std::add_const<value_type>::type>::type,
        reference>::type							const_reference;

      //typedef typename std::result_of<UnaryFunction(mln_creference(I))>::type	const_reference;

      struct pixel_type;
      struct const_pixel_type;


    private:
      struct val_fun_t {
        reference operator() (mln_reference(I) v) const {
          return m_fun(std::forward<mln_reference(I)>(v));
        }
        UnaryFunction m_fun;
      };

      struct const_val_fun_t {
        const_val_fun_t(const UnaryFunction& fun) : m_fun(fun)
        {
        }

        const_val_fun_t(const val_fun_t& other) : m_fun (other.m_fun)
        {
        }

        const_reference operator() (mln_reference(I) v) const {
          return m_fun(std::forward<mln_reference(I)>(v));
        }
        UnaryFunction m_fun;
      };

      struct pix_fun_t {
        pixel_type operator() (const mln_pixel(I)& px) const{
          return pixel_type(px, m_ima);
        }
        this_t* m_ima;
      };

      struct const_pix_fun_t {
        const_pix_fun_t(const this_t* ima)
        : m_ima(ima)
        {
        }

        const_pix_fun_t(const pix_fun_t& other)
        : m_ima (other.m_ima)
        {
        }

        const_pixel_type operator() (const mln_pixel(I)& px) const {
          return const_pixel_type(px, m_ima);
        }
        const this_t* m_ima;
      };


    public:
      // Ranges
      // \{
      typedef transformed_range<typename image_value_range<image_t>::type, val_fun_t>		value_range;
      typedef transformed_range<typename image_value_range<image_t>::type, const_val_fun_t>	const_value_range;
      typedef transformed_range<typename image_pixel_range<image_t>::type, pix_fun_t>		pixel_range;
      typedef transformed_range<typename image_pixel_range<image_t>::type, const_pix_fun_t>	const_pixel_range;
      // \}

      transformed_image() = default;

      transformed_image(I&& ima, const UnaryFunction& f)
        : transformed_image_helper<I> { std::forward<I>(ima) },
          m_fun(f)
      {
      }


      friend
      internal::initializer<mln_concrete(transformed_image),
                            typename internal::image_init_from<transformed_image>::type>
      imconcretize(const transformed_image& f)
      {
        using mln::imchvalue;
        return std::move(imchvalue<value_type>(f.m_ima));
      }

      template <typename V>
      friend
      internal::initializer<mln_ch_value(transformed_image, V),
                            typename internal::image_init_from<transformed_image>::type>
      imchvalue(const transformed_image& f)
      {
        using mln::imchvalue;
        return std::move(imchvalue<V>(f.m_ima));
      }


      const_value_range
      values() const
      {
        const_val_fun_t f {m_fun};
        return const_value_range(this->m_ima.values(), f);
      }

      value_range
      values()
      {
        val_fun_t f {m_fun};
        return value_range(this->m_ima.values(), f);
      }

      const_pixel_range
      pixels() const
      {
        const_pix_fun_t fun = { this };
        return const_pixel_range(this->m_ima.pixels(), fun);
      }

      pixel_range
      pixels()
      {
        pix_fun_t fun = { this };
        return pixel_range(this->m_ima.pixels(), fun);
      }



      template <typename dummy = reference>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      operator() (const mln_point(I)& p)
      {
        mln_precondition(this->domain().has(p));
        return m_fun(this->m_ima(p));
      }

      template <typename dummy = const_reference>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      operator() (const mln_point(I)& p) const
      {
        mln_precondition(this->domain().has(p));
        return m_fun(this->m_ima(p));
      }

      template <typename dummy = reference>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      at (const mln_point(I)& p)
      {
        return m_fun(this->m_ima.at(p));
      }

      template <typename dummy = const_reference>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      at (const mln_point(I)& p) const
      {
        return m_fun(this->m_ima.at(p));
      }

      template <typename dummy = pixel_type>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      pixel (const mln_point(I)& p)
      {
        mln_precondition(this->domain().has(p));
        return pixel_type(this->m_ima.pixel(p), this);
      }

      template <typename dummy = const_pixel_type>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      pixel (const mln_point(I)& p) const
      {
        mln_precondition(this->domain().has(p));
        return const_pixel_type(this->m_ima.pixel(p), this);
      }

      template <typename dummy = pixel_type>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      pixel_at (const mln_point(I)& p)
      {
        return pixel_type(this->m_ima.pixel_at(p), this);
      }

      template <typename dummy = const_pixel_type>
      typename std::enable_if< image_traits<this_t>::accessible::value, dummy >::type
      pixel_at (const mln_point(I)& p) const
      {
        return const_pixel_type(this->m_ima.pixel_at(p), this);
      }


      template <typename dummy = reference, typename image_type = image_t>
      typename std::enable_if< image_traits<this_t>::indexable::value, dummy >::type
      operator[] (typename image_type::size_type i)
      {
        return m_fun(this->m_ima[i]);
      }

      template <typename dummy = const_reference, typename image_type = image_t>
      typename std::enable_if< image_traits<this_t>::indexable::value, dummy >::type
      operator[] (typename image_type::size_type i) const
      {
        return m_fun(this->m_ima[i]);
      }

    private:
      UnaryFunction m_fun;
    };


    /*****************************/
    /**   Pixel definition      **/
    /*****************************/

    template <typename I, class UnaryFunction>
    struct transformed_image<I, UnaryFunction, false>::pixel_type
      : morpher_pixel_base<transformed_image<I, UnaryFunction, false>::pixel_type, mln_pixel(I)>
    {
      friend struct mln::morpher_core_access;
      friend struct const_pixel_type;
      typedef transformed_image<I, UnaryFunction, false>                           image_type;
      typedef typename transformed_image<I, UnaryFunction, false>::reference       reference;
      typedef typename transformed_image<I, UnaryFunction, false>::value_type      value_type;

      pixel_type(const mln_pixel(I)& px, image_type* ima)
	: m_pix(px), m_ima(ima)
      {
      }

      reference val() const
      {
	return m_ima->m_fun(m_pix.val());
      }

      image_type&	image() const
      {
	return *m_ima;
      }


    private:
      mln_pixel(I)		m_pix;
      image_type*		m_ima;
    };

    template <typename I, class UnaryFunction>
    struct transformed_image<I, UnaryFunction, false>::const_pixel_type
      : morpher_pixel_base<transformed_image<I, UnaryFunction, false>::const_pixel_type, mln_pixel(I)>
    {
      friend struct mln::morpher_core_access;
      typedef const transformed_image<I, UnaryFunction, false>                  image_type;
      typedef transformed_image<I, UnaryFunction, false>::value_type		value_type;
      typedef transformed_image<I, UnaryFunction, false>::const_reference	reference;

      const_pixel_type(const mln_pixel(I)& px, image_type* ima)
	: m_pix(px), m_ima(ima)
      {
      }

      // Interop
      const_pixel_type(const pixel_type& other)
      : m_pix(other.m_pix), m_ima(other.m_ima)
      {
      }


      reference val() const
      {
	return m_ima->m_fun(m_pix.val());
      }

      image_type&	image() const
      {
	return *m_ima;
      }


    private:
      mln_pixel(I)	   m_pix;
      image_type*	   m_ima;
    };

  }


  template <typename I, class UnaryFunction>
  inline
  transformed_image<I, UnaryFunction>
  imtransform(Image<I>&& ima, const UnaryFunction& f)
  {
    return transformed_image<I, UnaryFunction>(move_exact<I>(ima), f);
  }

  template <typename I, class UnaryFunction>
  transformed_image<const I&, UnaryFunction>
  imtransform(const Image<I>& ima, const UnaryFunction& f)
  {
    return transformed_image<const I&, UnaryFunction>(exact(ima), f);
  }

  template <typename I, class UnaryFunction>
  transformed_image<I&, UnaryFunction>
  imtransform(Image<I>& ima, const UnaryFunction& f)
  {
    return transformed_image<I&, UnaryFunction>(exact(ima), f);
  }

} // end of namespace mln

#endif //!MLN_CORE_IMAGE_MORPHERS_TRANSFORMED_IMAGE_HPP
