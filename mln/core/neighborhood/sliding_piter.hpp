#ifndef MLN_CORE_NEIGHBORHOOD_SLIDING_PITER_HPP
# define MLN_CORE_NEIGHBORHOOD_SLIDING_PITER_HPP

# include <mln/core/range/range.hpp>
# include <mln/core/iterator/iterator_base.hpp>

namespace mln
{

  /// \brief Define an iterator other a siteset centered on a site.
  ///
  /// Define an iterator other a siteset centered on a site.
  ///
  /// \p Point can be a pointer type, in that case, the current iterator will
  /// be binded to the point.
  /// \p Point can be an iterator, in that case, the current iterator will
  /// be binded to the iterator.
  ///
  /// Note that the siteset is not copied and its lifetime should be longer that
  /// the iterator's one.
  ///
  /// \tparam Point
  /// \tparam SiteSet
  template <class Point, class SiteSet, class Enable = void>
  struct sliding_piter;


  /******************************************/
  /****          Facade                  ****/
  /******************************************/


  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <class Point, class SiteSet, class Enable>
  struct sliding_piter :
    iterator_base< sliding_piter<Point, SiteSet>, Point, Point>
  {
    sliding_piter() = default;
    sliding_piter(const Point& p, const SiteSet& s)
      : m_p (p), m_it (rng::iter(s))
    {
    }

    void init() { m_it.init(); }
    void next() { m_it.next(); }
    bool finished() const { return m_it.finished(); }
    Point dereference() const { return m_p + *m_it; }


    void rebind(const Point& p) { m_p = p; }

  private:
    Point m_p;
    typename range_const_iterator<SiteSet>::type m_it;
  };


  template <class P, class SiteSet>
  struct sliding_piter<P, SiteSet,
		       typename std::enable_if<std::is_pointer<P>::value or
					       mln::is_a<P, Iterator>::value>::type
		       >
  : iterator_base< sliding_piter<P, SiteSet>,
		   typename range_value<SiteSet>::type,
		   typename range_value<SiteSet>::type >
  {
    sliding_piter(const P& p, const SiteSet& s)
      : m_p (p), m_it (rng::iter(s))
    {
    }

    void init() { m_it.init(); }
    void next() { m_it.next(); }
    bool finished() const { return m_it.finished(); }

    typename range_value<SiteSet>::type
    dereference() const { return *m_p + *m_it; }
    void rebind(const P& p) { m_p = &p; }

  private:
    typename std::conditional<std::is_pointer<P>::value,
			      const P, const P&>::type m_p;
    typename range_const_iterator<SiteSet>::type m_it;
  };



}

#endif // ! MLN_CORE_NEIGHBORHOOD_SLIDING_PITER_HPP
