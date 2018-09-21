#ifndef SHAPE_HPP
# define SHAPE_HPP

# include <type_traits>
# include <functional>
//# include <vector>
# include <unordered_set>
# include <ostream>
# include <mln/core/domain/box.hpp>

# include <boost/dynamic_bitset.hpp>
# include <boost/range/iterator_range.hpp>

/**
*
* \tparam LowerCompare Model of relation
* \tparam UpperCompare Model of relation
*
* A relation R must verify:
*
* \f[
*  a R b ^ b R c ⇒ a R c
*  a R b ^ a R c ⇒ a R k  V  a = k with k=relation_traits<R>::comm(b,c)
*
*
* \f]
*/
template <typename V, class LowerCompare, class UpperCompare>
struct shape;

template <typename V, class LowerCompare, class UpperCompare>
bool
is_shape_included(const shape<V, LowerCompare, UpperCompare>& a,
		  const shape<V, LowerCompare, UpperCompare>& b);

template <typename V, class LowerCompare, class UpperCompare>
bool
is_shape_equal(const shape<V, LowerCompare, UpperCompare>& a,
	       const shape<V, LowerCompare, UpperCompare>& b);



template <typename V, class LowerCompare, class UpperCompare>
struct shape
{
  enum { NONE = 0, LOWER = 1, UPPER = 2, LEQ = 4, GEQ = 8 } ;

private:
  static constexpr bool hasUL = not std::is_same<LowerCompare, UpperCompare>::value;

public:
  shape(V level, LowerCompare cmp, int sz, int ncols);

  template <class dummy = void>
  shape(V level, UpperCompare cmp, int sz, int ncols,
	typename std::enable_if<hasUL, dummy>::type* = NULL);

  shape(const shape& ) = delete;
  shape(shape&&);
  shape& operator= (const shape& ) = delete;
  shape& operator= (shape&&);

  bool operator== (const shape&) const;

  bool		islower() const;
  bool		isupper() const;
  bool		isleq() const;
  bool		isgeq() const;
  V		lower_level() const;
  V		upper_level() const;

  void		set_level(V level, LowerCompare cmp);

  template <typename dummy = void>
  void		set_level(V level, UpperCompare cmp,
			  typename std::enable_if<hasUL, dummy>::type* = NULL);

  void		init();
  void		add_point(const mln::point2d& p);


  unsigned		size() const;
  const mln::box2d&	bbox() const;


  // Update a shape with another shape which is equal.
  void		update_with(const shape& other) const;

  struct set_iterator;

  boost::iterator_range<set_iterator> pset() const;

private:
  friend bool is_shape_included<>(const shape<V, LowerCompare, UpperCompare>& a,
				const shape<V, LowerCompare, UpperCompare>& b);

  friend bool is_shape_equal<>(const shape<V, LowerCompare, UpperCompare>& a,
			       const shape<V, LowerCompare, UpperCompare>& b);


  mutable int				m_type;
  mutable V				m_level_1;
  mutable V				m_level_2;
  mutable LowerCompare			m_cmp_1;
  mutable UpperCompare			m_cmp_2;

  const size_t				m_size; ///< Number of points in the domain Ω
  const size_t			        m_nc;	///< Number of columns

  mln::box2d				m_bbox;   ///< BBox of the set
  boost::dynamic_bitset<>		m_set;    ///< BitSet of points
  boost::dynamic_bitset<>::size_type    m_first;  ///< First point in the set
  boost::dynamic_bitset<>::size_type	m_count;  ///< Number of points in the set

};

template <typename shape_t>
using shape_set = std::unordered_set<shape_t>;

template <typename V, class LowerCompare, class UpperCompare>
void
prettyprint_shape(const shape<V, LowerCompare, UpperCompare>& shp,
		  const mln::box2d& domain);

namespace std {
  /// \brief specialization of hash function for the unordered set of
  template <typename V, class LowerCompare, class UpperCompare>
  struct hash< shape<V, LowerCompare, UpperCompare> >;
}


/*************************************/
/***       Implementation          ***/
/*************************************/


template <typename V, class LowerCompare, class UpperCompare>
struct shape<V,LowerCompare,UpperCompare>::set_iterator :
  boost::iterator_facade<shape::set_iterator,
			 mln::point2d,
			 boost::forward_traversal_tag,
			 mln::point2d>
{
  set_iterator(const boost::dynamic_bitset<>* bset, int ncols,
	       boost::dynamic_bitset<>::size_type first)
    : m_bset(bset),
      m_nc(ncols),
      m_cur(first)
  {
  }


private:
  friend class boost::iterator_core_access;

  void increment()
  {
    m_cur = m_bset->find_next(m_cur);
  }

  bool equal(const set_iterator& other) const
  {
    return this->m_cur == other.m_cur;
  }

  mln::point2d dereference() const
  {
    return { (short)(m_cur / m_nc), (short)(m_cur % m_nc) };
  }


private:
  const boost::dynamic_bitset<>*	m_bset;
  int					m_nc;
  boost::dynamic_bitset<>::size_type	m_cur;
};



namespace impl
{
  inline
  bool
  is_box_inc(const mln::box2d& a, const mln::box2d& b)
  {
    return
      vecprod_isgreaterequal(a.pmin, b.pmin) and
      vecprod_islessequal(a.pmax, b.pmax);
  }

}

template <typename V, class LowerCompare, class UpperCompare>
bool
is_shape_included(const shape<V, LowerCompare, UpperCompare>& a,
		  const shape<V, LowerCompare, UpperCompare>& b)
{
  // Heuristique 1.
  // A c B => |A| < |B| et bbox(A) c bbox(B)
  if (not (a.size() <= b.size() and impl::is_box_inc(a.bbox(), b.bbox())))
    return false;

  // Heuristique 2.
  // A ∈ CC([u < a]), B ∈ CC([u < b]
  // a ≤ b ⇒ A ⋂ B = ∅ or A ⊂ B
  if (a.islower() and b.islower())
    {
      LowerCompare cmp;
      if (a.lower_level() == b.lower_level() or
	  cmp(a.lower_level(), b.lower_level()))
	return b.m_set.test(a.m_first);
    }
  if (a.isupper() and b.isupper() and
      !std::is_same<LowerCompare, UpperCompare>::value)
    {
      UpperCompare cmp;
      if (a.upper_level() == b.upper_level() or
	  cmp(a.upper_level(), b.upper_level()))
	return b.m_set.test(a.m_first);
    }

  // Heuristique 3.
  return a.m_set.is_subset_of(b.m_set);
}


template <typename V, class LowerCompare, class UpperCompare>
bool
is_shape_equal(const shape<V, LowerCompare, UpperCompare>& a,
	       const shape<V, LowerCompare, UpperCompare>& b)
{
  // Heuristique 1.
  // A = B => |A| = |B| et bbox(A) = bbox(B)
  if (not (a.size() == b.size() and a.bbox() == b.bbox()))
    return false;

  // Heuristique 2.
  // La < Lb or Lb < a => A ^ B = 0 or A c B or B c A
  // Heuristique 2.
  // La < Lb => A ^ B = 0 or A c B
  if (a.islower() and b.islower())
    {
      LowerCompare cmp;
      if (a.lower_level() == b.lower_level() or
	  cmp(a.lower_level(), b.lower_level()))
	return b.m_first == a.m_first;
    }
  if (a.isupper() and b.isupper() and
      !std::is_same<LowerCompare, UpperCompare>::value)
    {
      UpperCompare cmp;
      if (a.upper_level() == b.upper_level() or
	  cmp(a.upper_level(), b.upper_level()))
	return b.m_first == a.m_first;
    }

  // Heuristique 3.
  return b.m_set == a.m_set;

  //return std::equal(a.m_set.begin(), a.m_set.end(), b.m_set.begin());
}


template <typename V, class LowerCompare, class UpperCompare>
shape<V, LowerCompare, UpperCompare>::shape(V x, LowerCompare,
					    int sz, int ncols)
  : m_type(LOWER),
    m_level_1(x),
    m_size (sz),
    m_nc (ncols)
{
  init();
}

template <typename V, class LowerCompare, class UpperCompare>
template <class dummy>
shape<V, LowerCompare, UpperCompare>::shape(V x, UpperCompare,
					    int sz, int ncols,
					    typename std::enable_if<hasUL, dummy>::type*)
  : m_type(UPPER),
    m_level_2(x),
    m_size (sz),
    m_nc (ncols)
{
  init();
}

template <typename V, class LowerCompare, class UpperCompare>
shape<V, LowerCompare, UpperCompare>::shape(shape&& other)
  : m_type(other.m_type),
    m_level_1(other.m_level_1),
    m_level_2(other.m_level_2),
    m_cmp_1(other.m_cmp_1),
    m_cmp_2(other.m_cmp_2),
    m_size (other.m_size),
    m_nc (other.m_nc),
    m_bbox (other.m_bbox),
    m_first (other.m_first),
    m_count (other.m_count)
{
  m_set.swap(other.m_set);

  other.m_set.clear();
  other.m_type = NONE;
  other.m_first = m_set.npos;
  other.m_count = 0;
}

template <typename V, class LowerCompare, class UpperCompare>
shape<V, LowerCompare, UpperCompare>&
shape<V, LowerCompare, UpperCompare>::operator= (shape&& other)
{
  m_type = other.m_type;
  m_level_1 = other.m_level_1;
  m_level_2 = other.m_level_2;
  m_cmp_1 = other.m_cmp_1;
  m_cmp_2 = other.m_cmp_2;
  m_bbox = other.m_bbox;
  m_first = other.m_first;
  m_count = other.m_count;
  m_set.swap(other.m_set);

  other.m_set.clear();
  other.m_type = NONE;
  other.m_first = m_set.npos;
  other.m_count = 0;

  return *this;
}



template <typename V, class LowerCompare, class UpperCompare>
void
shape<V, LowerCompare, UpperCompare>::init()
{
  m_set.resize(m_size);
  m_set.reset();
  m_first = m_set.npos;
  m_count = 0;
  m_bbox.pmin = mln::value_traits<mln::point2d, mln::productorder_less<mln::point2d> >::sup();
  m_bbox.pmax = mln::value_traits<mln::point2d, mln::productorder_less<mln::point2d> >::inf();
}



template <typename V, class LowerCompare, class UpperCompare>
void
shape<V, LowerCompare, UpperCompare>::add_point(const mln::point2d& p)
{
  size_t x = p[0] * m_nc + p[1];

  if (m_count == 0) m_first = x;
  m_set.set(x);
  m_count++;

  m_bbox.pmin = mln::inf(m_bbox.pmin, p, mln::productorder_less<mln::point2d> ());
  m_bbox.pmax = mln::sup(m_bbox.pmax, p, mln::productorder_less<mln::point2d> ());
}

template <typename V, class LowerCompare, class UpperCompare>
void
shape<V, LowerCompare, UpperCompare>::set_level(V level, LowerCompare)
{
  m_type = LOWER;
  m_level_1 = level;
}


template <typename V, class LowerCompare, class UpperCompare>
template <typename dummy>
void
shape<V, LowerCompare, UpperCompare>::set_level(V level, UpperCompare,
						typename std::enable_if<hasUL, dummy>::type*)
{
  m_type = UPPER;
  m_level_2 = level;
}

template <typename V, class LowerCompare, class UpperCompare>
bool
shape<V, LowerCompare, UpperCompare>::operator== (const shape& other) const
{
  return is_shape_equal(*this, other);
}

template <typename V, class LowerCompare, class UpperCompare>
bool
shape<V, LowerCompare, UpperCompare>::isupper() const
{
  return (m_type & UPPER) != 0;
}

template <typename V, class LowerCompare, class UpperCompare>
bool
shape<V, LowerCompare, UpperCompare>::islower() const
{
  return (m_type & LOWER) != 0;
}

template <typename V, class LowerCompare, class UpperCompare>
bool
shape<V, LowerCompare, UpperCompare>::isleq() const
{
  return (m_type & LEQ) != 0;
}

template <typename V, class LowerCompare, class UpperCompare>
bool
shape<V, LowerCompare, UpperCompare>::isgeq() const
{
  return (m_type & GEQ) != 0;
}


template <typename V, class LowerCompare, class UpperCompare>
V
shape<V, LowerCompare, UpperCompare>::upper_level() const
{
  return m_level_2;
}

template <typename V, class LowerCompare, class UpperCompare>
V
shape<V, LowerCompare, UpperCompare>::lower_level() const
{
  return m_level_1;
}

template <typename V, class LowerCompare, class UpperCompare>
unsigned
shape<V, LowerCompare, UpperCompare>::size() const
{
  return m_count;
}

template <typename V, class LowerCompare, class UpperCompare>
const mln::box2d&
shape<V, LowerCompare, UpperCompare>::bbox() const
{
  return m_bbox;
}

template <typename V, class LowerCompare, class UpperCompare>
boost::iterator_range<typename shape<V, LowerCompare, UpperCompare>::set_iterator>
shape<V, LowerCompare, UpperCompare>::pset() const
{
  return boost::make_iterator_range(set_iterator(&m_set, m_nc, m_first),
				    set_iterator(&m_set, m_nc, m_set.npos));
}



template <typename V, class LowerCompare, class UpperCompare>
void
shape<V, LowerCompare, UpperCompare>::update_with(const shape& other) const
{
  mln_precondition(this->m_type != NONE and other.m_type != NONE);

  // Same shape type [u ≤ a] ∧ [u ≤ b] => [u ≤ inf(a,b)]
  if (other.islower()) {
    if (!this->islower())
      m_level_1 = other.m_level_1;
    else {
      m_level_1 = m_cmp_1.inf(other.m_level_1, this->m_level_1);
      m_type |= LEQ;
    }
  }

  // Same shape type [u ≥ a] ^ [u ≥ b] => [u ≥ sup(a,b)]
  if (other.isupper()) {
    if (!this->isupper())
      m_level_2 = other.m_level_2;
    else {
      m_level_2 = m_cmp_2.inf(other.m_level_2, this->m_level_2);
      m_type |= GEQ;
    }
  }
  m_type |= other.m_type;
}

template <typename V, class LowerCompare, class UpperCompare>
std::ostream&
operator<< (std::ostream& os, const shape<V, LowerCompare, UpperCompare>& shp)
{
  os << "==Info== Size:" << shp.size();
  if (shp.islower()) {
    os << "\tfrom: [u" << LowerCompare::str
       << (shp.isleq() ? '=' : ' ');
    format(os, shp.lower_level()) << "]";
  }
  if (shp.isupper()) {
    os << "\tfrom: [u" << UpperCompare::str
       << (shp.isgeq() ? '=' : ' ');
    format(os, shp.upper_level()) << "]";
  }
  return os;
}



template <typename V, class LowerCompare, class UpperCompare>
void
prettyprint_shape(const shape<V, LowerCompare, UpperCompare>& shp,
		  const mln::box2d& domain)
{
  using namespace mln;
  using mln::io::format;

  image2d<bool> ima;
  ima.resize(domain, 0, false);

  for (auto x: shp.pset())
    ima(x) = true;

  // if (not std::is_same<Ktag, K0_tag>::value)
  //   for (auto x: shp.nonprimaries)
  // 	ima(x) = true;

  std::cout << "==Info== Size:" << shp.size() <<
    "(0x" << std::hex << shp.size() << ")" << std::dec
	    << std::endl;

  if (shp.islower()) {
    std::cout << "\tfrom: [u" << LowerCompare::str
	      << (shp.isleq() ? '=' : ' ');
    format(std::cout, shp.lower_level()) << "]" << std::endl;
  }

  if (shp.isupper()) {
    std::cout << "\tfrom: [u" << UpperCompare::str
	      << (shp.isgeq() ? '=' : ' ');
    format(std::cout, shp.upper_level()) << "]" << std::endl;
  }

  std::cout << "hash: " << std::hex
	    << "0x" << std::hash< shape<V, LowerCompare, UpperCompare> >() (shp)
	    << std::dec << std::endl;

  io::imprint(ima);
}


namespace std
{


  template <typename V, class LowerCompare, class UpperCompare>
  struct hash< shape<V, LowerCompare, UpperCompare> >
  {
    typedef shape<V, LowerCompare, UpperCompare> shape_t;

    size_t
    operator() (const shape_t& x) const
    {
      size_t h = ((x.size() & 0x00FF) << 24) |
	((x.bbox().pmin[0] & 0x00FF) << 16) |
	((x.bbox().pmin[1] & 0x00FF) << 8) |
	((x.bbox().pmax[0] + x.bbox().pmax[1]) & 0x00FF);
      return h;
    }
  };

}

#endif // ! SHAPE_HPP
