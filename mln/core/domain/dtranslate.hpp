#ifndef MLN_CORE_DOMAINS_DTRANSLATE_HPP
# define MLN_CORE_DOMAINS_DTRANSLATE_HPP

# include <mln/core/range/adaptor/transformed.hpp>

namespace mln {

  namespace internal
  {
    template <typename Domain>
    struct translate_fun
    {
      typedef typename Domain::value_type P;
      typedef std::binder1st< std::plus<P> > type;
    };

  };


  template <typename Domain>
  using translated_domain = adaptors::transformed_range<const Domain&, typename internal::translate_fun<Domain>::type>;

  template <typename Domain>
  inline
  translated_domain<Domain>
  dtranslate(const Domain& domain, typename Domain::value_type p)
  {
    typedef typename Domain::value_type P;
    return adaptors::transform(domain, std::bind1st(std::plus<P> (), p));
  }


} // end of namespace mln

#endif //!MLN_CORE_DOMAINS_DTRANSLATE_HPP
