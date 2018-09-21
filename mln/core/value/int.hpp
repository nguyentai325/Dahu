#ifndef MLN_CORE_VALUE_INT_HPP
# define MLN_CORE_VALUE_INT_HPP

# include <cstdint>
# include <boost/integer.hpp>

# include <mln/core/value/value_traits.hpp>
# include <mln/core/value/index.hpp>


namespace mln
{

  typedef uint8_t       uint8;
  typedef uint16_t      uint16;
  typedef uint32_t      uint32;
  typedef uint64_t      uint64;

  typedef int8_t        int8;
  typedef int16_t       int16;
  typedef int32_t       int32;
  typedef int64_t       int64;

  template <unsigned nbits>
  struct UInt;

  template <unsigned nbits>
  struct Int;

} // end of namespace mln

  /******************************************/
  /****              traits              ****/
  /******************************************/

namespace std
{
  template <unsigned nbits>
  struct is_arithmetic< mln::UInt<nbits> > : std::true_type
  {
  };

  template <unsigned nbits>
  struct is_arithmetic< mln::Int<nbits> > : std::true_type
  {
  };

  template <unsigned nbits>
  struct is_signed< mln::Int<nbits> > : std::true_type
  {
  };

  template <unsigned nbits>
  struct is_unsigned< mln::UInt<nbits> > : std::true_type
  {
  };


  template <unsigned nbits>
  struct make_signed< mln::UInt<nbits> >
  {
    typedef mln::Int<nbits> type;
  };

  template <unsigned nbits>
  struct make_signed< mln::Int<nbits> >
  {
    typedef mln::Int<nbits> type;
  };


  template <unsigned nbits>
  struct make_unsigned< mln::UInt<nbits> >
  {
    typedef mln::UInt<nbits> type;
  };

  template <unsigned nbits>
  struct make_unsigned< mln::Int<nbits> >
  {
    typedef mln::UInt<nbits> type;
  };

}

namespace mln
{

  // Note that typename std::enable_if<std::is_arithmetic< UInt<nbits> >::value>::type
  // instead of just void is required to avoid ambiguous specialization
  template <unsigned nbits>
  struct value_traits<UInt<nbits>, std::less< UInt<nbits> >,
		      typename std::enable_if<std::is_arithmetic< UInt<nbits> >::value>::type >
  {
    static constexpr unsigned ndim = 1;
    static constexpr unsigned quant = nbits;
    static constexpr UInt<nbits> min() { return 0; }
    static constexpr UInt<nbits> max() { return (1ul << nbits) - 1; }
    static constexpr UInt<nbits> inf() { return min(); }
    static constexpr UInt<nbits> sup() { return max(); }
  };

  template <unsigned nbits>
  struct value_traits<Int<nbits>, std::less< Int<nbits> >,
		      typename std::enable_if<std::is_arithmetic< Int<nbits> >::value>::type >
  {
    static constexpr unsigned ndim = 1;
    static constexpr unsigned quant = nbits;
    static constexpr Int<nbits> min() { return -(1 << (nbits-1)); }
    static constexpr Int<nbits> max() { return (1ul << (nbits-1)) - 1; }
    static constexpr Int<nbits> inf() { return min(); }
    static constexpr Int<nbits> sup() { return max(); }
  };

} // end of namespace mln

  /******************************************/
  /****     Indexers specialization     *****/
  /******************************************/

namespace mln {

  // Forward declaration
  template <typename V, typename StrictWeakOrdering, typename Enable>
  struct indexer;

  template <unsigned nbits>
  struct indexer<UInt<nbits>, std::less< UInt<nbits> >, void>
  {
    typedef Index<UInt<nbits>, std::less< UInt<nbits> > > index_type;
    typedef UInt<nbits> value_type;
    static const bool inversible = true;

    index_type operator() (value_type x) const { return index_type(x); }
    value_type inv(index_type i)         const { return i; }
  };

  template <unsigned nbits>
  struct indexer<UInt<nbits>, std::greater< UInt<nbits> >, void>
  {
    typedef Index<UInt<nbits>, std::greater< UInt<nbits> > > index_type;
    typedef UInt<nbits> value_type;
    static const bool inversible = true;

    index_type operator() (value_type x) const { return index_type(x); }
    value_type inv(index_type i)         const { return i; }
  };

  template <unsigned nbits>
  struct indexer<Int<nbits>, std::less< Int<nbits> >, void>
  {
    typedef Index<Int<nbits>, std::less< Int<nbits> > > index_type;
    typedef Int<nbits> value_type;
    static const bool inversible = true;

    index_type operator() (value_type x) const { return index_type(x); }
    value_type inv(index_type i)         const { return i; }
  };

  template <unsigned nbits>
  struct indexer<Int<nbits>, std::greater< Int<nbits> >, void>
  {
    typedef Index<Int<nbits>, std::less< Int<nbits> > > index_type;
    typedef Int<nbits> value_type;
    static const bool inversible = true;

    index_type operator() (value_type x) const { return index_type(x); }
    value_type inv(index_type i)         const { return i; }
  };



  /******************************************/
  /****          implementation          ****/
  /******************************************/

  template <unsigned nbits>
  struct UInt
  {
  private:
    typedef typename boost::uint_t<nbits>::least enc;

  public:
    UInt() = default;
    UInt(enc x) : m_x(x) {}
    operator enc& ()          { return m_x; }
//operator enc  () const     { return m_x; }
    operator unsigned () const { return m_x; }

  private:
    enc m_x;
  };

  template <unsigned nbits>
  struct Int
  {
  private:
    typedef typename boost::int_t<nbits>::least enc;

  public:
    Int() = default;
    Int(enc x) : m_x(x) {}
    operator enc& ()           { return m_x; }
//    operator enc  () const     { return m_x; }
    operator int () const { return m_x; }

  private:
    enc m_x;
  };


} // end of namespace mln

// Specialization of standards traits with Int and UInt

namespace std
{
  // Forward declaration
  template <typename... types_>
  struct common_type;


  template <unsigned n1, unsigned n2>
  struct common_type< mln::UInt<n1>, mln::UInt<n2> >
  {
    typedef mln::UInt<(n1 > n2 ? n1 : n2)> type;
  };

  template <unsigned n, typename V>
  struct common_type< mln::UInt<n>, V >
    : common_type< typename boost::uint_t<n>::least, V >
  {
  };

  template <unsigned n, typename V>
  struct common_type< V, mln::UInt<n> >
    : common_type< V, typename boost::uint_t<n>::least >
  {
  };


  template <unsigned n1, unsigned n2>
  struct common_type< mln::Int<n1>, mln::UInt<n2> >
  {
    typedef mln::Int<(n1 > n2 ? n1 : n2)> type;
  };

  template <unsigned n, typename V>
  struct common_type< mln::Int<n>, V >
    : common_type< typename boost::int_t<n>::least, V >
  {
  };

  template <unsigned n, typename V>
  struct common_type< V, mln::Int<n> >
    : common_type< V, typename boost::int_t<n>::least >
  {
  };


}

#endif //!MLN_CORE_VALUE_INT_HPP
