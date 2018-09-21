#ifndef MLN_CORE_RANGE_ZIP_HPP
# define MLN_CORE_RANGE_ZIP_HPP

# include <mln/core/range/range.hpp>
# include <mln/core/iterator/zip_iterator.hpp>
# include <mln/core/internal/tuple_utility.hpp>

namespace mln
{

  template <class RangeTuple>
  struct zip_range;


  /******************************/
  /**** Implementation         **/
  /******************************/

  namespace internal
  {

    struct meta_get_iterator
    {
      template <typename R>
      struct apply {
	typedef typename range_iterator<typename std::remove_reference<R>::type>::type type;
      };

      template <typename R>
      typename range_iterator<R>::type
      operator() (R& range) const
      {
	return rng::iter(range);
      }
    };

    struct meta_get_const_iterator
    {
      template <typename R>
      struct apply {
	typedef typename range_const_iterator<typename std::remove_reference<R>::type>::type type;
      };

      template <typename R>
      typename range_const_iterator<R>::type
      operator() (const R& range) const
      {
	return rng::iter(range);
      }
    };

  }

  template <class... TRanges>
  struct zip_range< std::tuple<TRanges...> >
  {
  private:
    template <typename>
    friend struct zip_range;

    typedef std::tuple<TRanges...> RangeTuple;

    typedef std::tuple<typename range_iterator<
			 typename std::remove_reference<TRanges>::type
			 >::type...>		iterator_tuple;
    typedef std::tuple<typename range_const_iterator<
			 typename std::remove_reference<TRanges>::type
			 >::type...>		const_iterator_tuple;

  public:
    typedef zip_iterator<iterator_tuple>		iterator;
    typedef zip_iterator<const_iterator_tuple>		const_iterator;
    typedef typename iterator::value_type		value_type;
    typedef typename iterator::reference		reference;

    zip_range(const std::tuple<TRanges...>& rng)
    : m_rng (rng)
    {
    }

    zip_range(const zip_range&) = default;

    //interop
    template <typename R2>
    zip_range(const zip_range<R2>& other,
	      typename std::enable_if<std::is_convertible<R2, RangeTuple>::value>::type* = NULL)
      : m_rng(other.m_rng)
    {
    }

    const_iterator iter() const
    {
      auto t = internal::tuple_transform(m_rng, internal::meta_get_const_iterator());
      return const_iterator(t);
    }

    iterator iter()
    {
      auto t = internal::tuple_transform(m_rng, internal::meta_get_iterator());
      return iterator(t);
    }

    std::size_t size() const
    {
      return m_rng.size();
    }

  private:
    RangeTuple		m_rng;
  };



  namespace rng
  {

    template <class... TRanges>
    zip_range< std::tuple<TRanges...> >
    zip(TRanges&&... ranges)
    {
      typedef std::tuple<TRanges...> T;
      return zip_range<T>(std::forward_as_tuple<TRanges&&...>(ranges...));
    }

  } // end of namespace mln::rng

} // end of namespace mln

#endif // ! MLN_CORE_RANGE_ZIP_HPP
