#ifndef MLN_CORE_RANGE_RBEGIN_HPP
# define MLN_CORE_RANGE_RBEGIN_HPP

namespace mln {

  template< class C >
  auto rbegin( C& c ) -> decltype(c.rbegin())
  {
    return c.rbegin();
  }

} // end of namespace mln

#endif //!MLN_CORE_RANGE_RBEGIN_HPP
