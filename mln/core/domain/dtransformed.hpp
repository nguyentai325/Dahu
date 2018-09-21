#ifndef MLN_CORE_DOMAINS_DTRANSFORMED_HPP
# define MLN_CORE_DOMAINS_DTRANSFORMED_HPP

# include <mln/core/range.hpp>

/// \file
/// \brief Lazy transformation of a domain through a function

namespace mln {

  // FWD
  namespace internal {
    template <typename Domain, typename Function, typename domain_category>
    struct internal::domain_transform;
  }

  /**
   * \brief Lazy application of the function to domain
   *
   * \param domain The domain to be mapped.
   * \param f The function to map.
   *
   * \post The following sementic is hold:
   * \code
   * for (x,y: zip(domain,  dtransform(domain, f))
   *    mln_invariant(y == f(x))
   * \endcode
   **/
  template <typename Domain, typename Function>
  internal::domain_transform<Domain, Function, typename range_category<Domain>::type>
  dtransform(const Domain& domain, Function f);


  //******************
  // Implementation  *
  //******************


  namespace internal
  {

    template <typename Domain, typename Function>
    struct domain_transform<Domain, Function, single_pass_range_tag>
    {
      typedef single_pass_range_tag                             category;
      typedef boost::range_value<Domain>::type                  value_type;
      typedef boost::range_iterator<Domain>::type               base_iterator_;
      typedef boost::transform_iterator<Function, base_iterator_> iterator;
      typedef iterator                                          const_iterator;

      domain_transform(const Domain& d, Function f):
        domain_ (d), f_ (f)
      {}

      iterator begin() const { return iterator(domain.begin(), f_); }
      iterator end() const { return iterator(domain.end(), f_); }

    protected:
      D domain_;
      Function f_;
    };

    template <typename Domain, typename Function>
    struct domain_transform<Domain, Function, forward_range_tag>
      : domain_transform<Domain, Function, single_pass_range_tag>
    {
      typedef domain_transform<Domain, Function, single_pass_range_tag> base;
      domain_transform(const Domain& d, Function f): base(d, f) {}
    }

    template <typename Domain, typename Function>
    struct domain_transform<Domain, Function, bidirectional_tag>
      : domain_transform<Domain, Function, forward_range_tag>
    {
      typedef boost::transform_iterator<Function, Domain>       iterator;
      typedef iterator                                          const_iterator;

      typedef domain_transform<Domain, Function, forward_range_tag> base;
      domain_transform(const Domain& d, Function f): base(d, f) {}

      iterator begin() const { return iterator(domain.begin(), f_); }
      iterator end() const { return iterator(domain.end(), f_); }
    }






  }


} // end of namespace mln


#endif //!MLN_CORE_DOMAINS_DTRANSFORMED_HPP
