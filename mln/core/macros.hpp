#ifndef MLN_CORE_MACROS_HPP
# define MLN_CORE_MACROS_HPP

// # include <boost/mpl/if.hpp>
// # include <boost/foreach.hpp>
// # include <utility>

// # define mln_site(I) typename I::site_type
// # define mln_point(I) typename I::point_type
// # define mln_value(I) typename I::value_type

// # define mln_piter(I) typename I::domain_type::iterator
// # define mln_viter(I) typename I::iterator
// # define mln_pixter(I) typename boost::mpl::if_< boost::is_const<I>, I::const_pixel_iterator, I::pixel_iterator >::type
// # define mln_pixter_(I) I::pixel_iterator

// # define mln_neighb(nbhor,I) typename typeof(nbhor)::nbh<I>::type
// # define mln_neighb_(nbhor,I) typeof(nbhor)::nbh<I>::type

// # define forallp(x, ima) BOOST_FOREACH(auto x, ima.domain())
// # define forallv(x, ima) BOOST_FOREACH(auto x, ima.values())

// # define forall(x, ima) BOOST_FOREACH(auto x, ima.pixels())

#endif /* !MLN_CORE_MACROS_HPP */
