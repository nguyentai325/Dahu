#ifndef MLN_CORE_ALGORITHM_SORT_INDEXES_HPP
# define MLN_CORE_ALGORITHM_SORT_INDEXES_HPP

# include <mln/core/value/value_traits.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/value/indexer.hpp>
# include <type_traits>
# include <vector>
# include <algorithm>

// FIXME: Speed up HQ version
// FIXME: Radix sort bugs after 32

namespace mln
{

  template <typename I, typename BinaryFunction = std::less<mln_value(I)> >
  std::vector<typename I::size_type>
  sort_indexes(const Image<I>& input,  BinaryFunction cmp = BinaryFunction ());


  template <typename I, typename OutputIterator, typename BinaryFunction = std::less<mln_value(I)> >
  void
  sort_indexes_it(const Image<I>& input, OutputIterator out, BinaryFunction cmp = BinaryFunction ());


  /****************/
  /* Implem       */
  /****************/

  namespace impl
  {
    struct use_counting_sort_tag  {};
    struct use_radix_sort_tag	  {};
    struct use_std_sort_tag	  {};


    template <typename I, typename OutputIterator, typename StrictWeakOrdering,
	      typename Indexer = indexer<mln_value(I), StrictWeakOrdering> >
    void
    sort_indexes(const I& input, OutputIterator v, StrictWeakOrdering, use_counting_sort_tag)
    {
      typedef typename Indexer::index_type index_t;
      Indexer f;

      static constexpr std::size_t nvalues = 1 << value_traits<index_t>::quant;
      unsigned h[nvalues] = {0,};
      {
	mln_viter(v, input);
	mln_forall(v)
	  ++h[f(*v)];

	unsigned count = 0;
        index_t i = value_traits<index_t>::min();
	do {
          unsigned tmp = h[i];
          h[i] = count;
          count += tmp;
        } while (i++ < value_traits<index_t>::max());
        assert(count == input.domain().size());
      }

      {
	mln_pixter(px, input);
	mln_forall(px)
	  v[ h[f(px->val())]++ ] = px->index();
      }
    }


    template <typename I, typename StrictWeakOrdering,
	      typename size_type = typename I::size_type,
	      typename Indexer = indexer<mln_value(I), StrictWeakOrdering> >
    void
    sort_indexes(const I& input, size_type* v, StrictWeakOrdering, use_radix_sort_tag)
    {
      typedef typename Indexer::index_type index_t;
      Indexer f;

      std::size_t n = input.domain().size();

      static constexpr std::size_t nvalues = 1 << 16;
      static constexpr std::size_t nvalues2 = 1 << (value_traits<index_t>::quant - 16);
      unsigned h[nvalues] = {0,};
      unsigned h2[nvalues2] = {0, };

      size_type* buffer;
      std::tie(buffer, std::ignore) = std::get_temporary_buffer<size_type>(n);

      // Last digit first
      {
	mln_viter(v, input);
	mln_forall(v) {
	  unsigned x = f(*v);
	  ++h[x & 0x0000FFFF];
	  ++h2[(x >> 16) & 0x0000FFFF];
	}

	{
	  unsigned count = 0;
	  index_t i( (unsigned) (value_traits<index_t>::min()) & 0x0000FFFF);
	  index_t j( (unsigned) (value_traits<index_t>::max()) & 0x0000FFFF);
	  do {
	    unsigned tmp = h[i];
	    h[i] = count;
	    count += tmp;
	  } while (i++ < j);
	}

	{
	  unsigned count = 0;
	  index_t i( ((unsigned) (value_traits<index_t>::min()) >> 16) & 0x0000FFFF);
	  index_t j( ((unsigned) (value_traits<index_t>::max()) >> 16) & 0x0000FFFF);
	  do {
	    unsigned tmp = h2[i];
	    h2[i] = count;
	    count += tmp;
	  } while (i++ < j);
	}


	mln_pixter(px, input);
	mln_forall(px) {
	  unsigned x = f(px->val());
	  x &= 0x0000FFFF;
	  buffer[ h[x]++ ] = px->index();
	}
      }

      // Next base 2**16 digit
      {
	unsigned shft = 16;
	while (true)
	  {
	    //std::fill(h, h + nvalues, 0);
	    // for (std::size_t i = 0; i < n; ++i) {
	    //   std::size_t p = buffer[i];
	    //   unsigned x = f(input[p]);
	    //   x = (x >> shft) & 0x0000FFFF;
	    //   ++h[x];
	    // }

	    // unsigned count = 0;
	    // index_t i( ((unsigned)(value_traits<index_t>::min()) >> shft) & 0x0000FFFF);
	    // index_t j( ((unsigned)(value_traits<index_t>::max()) >> shft) & 0x0000FFFF);
	    // do {
	    //   unsigned tmp = h[i];
	    //   h[i] = count;
	    //   count += tmp;
	    // } while (i++ < j);
	    // assert(count == n);

	    if (shft + 16 >= value_traits<index_t>::quant)
	      {
		for (std::size_t i = 0; i < n; ++i) {
		  std::size_t p = buffer[i];
		  unsigned x = f(input[p]);
		  x = (x >> shft) & 0x0000FFFF;
		  v[ h2[x]++ ] = p;
		}
		break;
	      }
	    else
	      {
		abort();
	      }
	  }
      }
      std::return_temporary_buffer(buffer);
    }


    template <typename I, typename OutputIterator, typename StrictWeakOrdering,
	      typename Indexer = indexer<mln_value(I), StrictWeakOrdering> >
    void
    sort_indexes(const I& input, OutputIterator v, StrictWeakOrdering cmp, use_std_sort_tag)
    {
      typedef typename I::size_type size_type;

      std::size_t i = 0;
      mln_pixter(px, input);
      mln_forall(px)
	v[i++] = px->index();

      std::sort(v, v + input.domain().size(), [&input, cmp](size_type x, size_type y) { return cmp(input[x], input[y]); });
    }


  } // end of namespace mln::impl


  template <typename I, typename BinaryFunction>
  std::vector<typename I::size_type>
  sort_indexes(const Image<I>& input, BinaryFunction cmp)
  {
    static_assert(std::is_same<typename image_category<I>::type, raw_image_tag>::value,
		  "Image must model the Raw Image Concept");
    typedef typename
      std::conditional< (value_traits<mln_value(I)>::quant <= 18), impl::use_counting_sort_tag, impl::use_radix_sort_tag >::type dispatch_tag;

    std::vector<typename I::size_type> v;
    v.resize(exact(input).domain().size());

    impl::sort_indexes(exact(input), &v[0], cmp, dispatch_tag ());
    return v;
  }

  template <typename I, typename OutputIterator, typename BinaryFunction>
  void
  sort_indexes_it(const Image<I>& input, OutputIterator out, BinaryFunction cmp)
  {
    static_assert(std::is_same<typename image_category<I>::type, raw_image_tag>::value,
		  "Image must model the Raw Image Concept");

    typedef typename
      std::conditional< (value_traits<mln_value(I)>::quant <= 18), impl::use_counting_sort_tag, impl::use_radix_sort_tag >::type dispatch_tag;

    impl::sort_indexes(exact(input), out, cmp, dispatch_tag ());
  }


} // end of namespace mln


#endif // !MLN_CORE_ALGORITHM_SORT_INDEXES_HPP
