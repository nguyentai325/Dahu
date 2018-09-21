#ifndef MLN_CORE_RANGE_REND_HPP
# define MLN_CORE_RANGE_REND_HPP

namespace mln {

  template< class C >
  auto rend( C& c ) -> decltype(c.rend())
  {
    return c.rend();
  }

} // end of namespace mln

#endif //!MLN_CORE_RANGE_REND_HPP
