// Copyright (C) 2009 EPITA Research and Development Laboratory (LRDE)
//
// This file is part of Olena.
//
// Olena is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, version 2 of the License.
//
// Olena is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Olena.  If not, see <http://www.gnu.org/licenses/>.
//
// As a special exception, you may use this file as part of a free
// software project without restriction.  Specifically, if other files
// instantiate templates or use macros or inline functions from this
// file, or you compile this file and link it with other files to produce
// an executable, this file does not by itself cause the resulting
// executable to be covered by the GNU General Public License.  This
// exception does not however invalidate any other reasons why the
// executable file might be covered by the GNU General Public License.

#ifndef MLN_MOPRHO_PSET_PSET_HPP
# define MLN_MOPRHO_PSET_PSET_HPP

# ifndef PSET_SWITCH_HQUEUE_NBITS
#  define PSET_SWITCH_HQUEUE_NBITS 18
# endif

# include <mln/core/value/value_traits.hpp>
# include <mln/core/value/indexer.hpp>
# include <set>
# include <algorithm>

namespace mln
{

  namespace morpho
  {

    //bool Enable = has_indexer<typename I::value_type, Compare>::value  >
    template <typename I,
	      typename Compare = std::less<typename I::value_type>,
	      bool Enable = (value_traits<mln_value(I)>::quant <= PSET_SWITCH_HQUEUE_NBITS and
			     has_indexer<mln_value(I), Compare>::value) >
    struct pset;


    template <typename I, typename Compare>
    struct pset<I, Compare, false>
    {
      typedef typename I::value_type	key_type;
      typedef typename I::size_type	value_type;

      struct cmp_t
      {
	bool operator () (const value_type& a, const value_type& b)
	{
	  return m_cmp(m_ima[a], m_ima[b]);
	}

	bool cmpkey (const key_type& x, const value_type& b)
	{
	  return m_cmp(x, m_ima[b]);
	}

	const I& m_ima;
	Compare  m_cmp;
      };

      pset(const I& ima, const Compare& cmp = Compare() );

      void insert(const value_type& v);
      bool empty() const;
      bool has_previous(const value_type& v) const;
      bool has_next(const value_type& v)	  const;

      value_type find_next(const value_type& v) const;
      value_type find_previous(const value_type& v) const;
      value_type pop_next(const value_type& v);
      value_type pop_previous(const value_type& v);

    private:
      std::multiset<std::size_t, cmp_t>		m_set;
    };


    template <typename I, typename Compare>
    struct pset<I, Compare, true>
    {
      typedef mln_value(I)		key_type;
      typedef typename I::size_type	value_type;

      pset(const I& ima, const Compare& cmp);

      void insert(const value_type& v);
      bool empty() const;
      bool has_previous(const value_type& v) const;
      bool has_next(const value_type& v)	  const;

      value_type find_next(const value_type& v) const;
      value_type find_previous(const value_type& v) const;
      value_type pop_next(const value_type& v);
      value_type pop_previous(const value_type& v);

    private:
      typedef typename I::size_type			size_type;
      typedef typename indexer<mln_value(I), Compare>::index_type  index_type;
      static constexpr size_type npos = -1;
      static constexpr std::size_t nvalues = indexer<mln_value(I), Compare>::nvalues;

      const I&		m_ima;
      Compare		m_cmp;
      indexer<mln_value(I), Compare> h;

      mln_ch_value(I, size_type) m_next;
      size_type			 m_hq[nvalues];
      unsigned			 m_size;
    };


    /*************************/
    /**  Implementation     **/
    /*************************/



    template <typename I, typename Compare>
    pset<I, Compare, false>::pset(const I& ima, const Compare& cmp)
      : m_set(cmp_t {ima, cmp} )
    {
    }

    template <typename I, typename Compare>
    void
    pset<I, Compare, false>::insert(const value_type& v)
    {
      m_set.insert(v);
    }

    template <typename I, typename Compare>
    bool
    pset<I, Compare, false>::empty() const
    {
      return m_set.empty();
    }

    template <typename I, typename Compare>
    bool
    pset<I, Compare, false>::has_previous(const value_type& p) const
    {
      mln_precondition(!empty());
      return m_set.key_comp()(*(m_set.begin()), p);
    }

    template <typename I, typename Compare>
    bool
    pset<I, Compare, false>::has_next(const value_type& p) const
    {
      mln_precondition(!empty());
      return !m_set.key_comp()(*(m_set.rbegin()), p);
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset<I, Compare, false>::find_previous(const value_type& p) const
    {
      mln_precondition(has_previous(p));
      return *m_set.lower_bound(p);
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset<I, Compare, false>::find_next(const value_type& p) const
    {
      mln_precondition(has_next(p));
      return *m_set.upper_bound(p);
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset<I, Compare, false>::pop_previous(const value_type& p)
    {
      mln_precondition(has_previous(p));
      auto x = m_set.upper_bound(p);
      --x;
      mln_assertion(x != m_set.end());
      value_type v = *x;
      m_set.erase(x);
      return v;
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset<I, Compare, false>::pop_next(const value_type& p)
    {
      mln_precondition(has_next(p));
      auto x = m_set.lower_bound(p);
      mln_assertion(x != m_set.end());
      value_type v = *x;
      m_set.erase(x);
      return v;
    }

    template <typename I, typename Compare>
    inline
    pset<I, Compare, true>::pset(const I& ima, const Compare& cmp)
      : m_ima(ima),
	m_cmp(cmp),
	m_size (0)
    {
      resize(m_next, m_ima).init((size_type)npos);
      std::fill(m_hq, m_hq + nvalues, (size_type)npos);
    }

    template <typename I, typename Compare>
    inline
    void
    pset<I, Compare, true>::insert(const value_type& v)
    {
      auto k = h(m_ima[v]);
      m_next[v] = m_hq[k];
      m_hq[k] = v;
      ++m_size;
    }

    template <typename I, typename Compare>
    inline
    bool
    pset<I, Compare, true>::empty() const
    {
      return m_size == 0;
    }


    template <typename I, typename Compare>
    inline
    bool
    pset<I, Compare, true>::has_previous(const value_type& v) const
    {
      for (auto k = h(m_ima[v]); k > 0;)
	if (m_hq[--k] != npos)
	  return true;
      return false;
    }

    template <typename I, typename Compare>
    inline
    bool
    pset<I, Compare, true>::has_next(const value_type& v) const
    {
      auto k = h(m_ima[v]);
      for (; k < value_traits<index_type>::max(); ++k)
	if (m_hq[k] != npos)
	  return true;
      return (m_hq[k] != npos);
    }

    template <typename I, typename Compare>
    inline
    typename pset<I, Compare, true>::value_type
    pset<I, Compare, true>::find_previous(const value_type& v) const
    {
      mln_precondition(!empty());
      auto k = h(m_ima[v]);
      while (m_hq[--k] == npos)
	;

      mln_postcondition(k < h(m_ima[v]));
      return m_hq[k];
    }


    template <typename I, typename Compare>
    inline
    typename pset<I, Compare, true>::value_type
    pset<I, Compare, true>::find_next(const value_type& v) const
    {
      mln_precondition(!empty());
      auto k = h(m_ima[v]);
      while (m_hq[k] == npos)
	++k;

      mln_postcondition(not (k < h(m_ima[v])));
      return m_hq[k];
    }

    template <typename I, typename Compare>
    inline
    typename pset<I, Compare, true>::value_type
    pset<I, Compare, true>::pop_previous(const value_type& v)
    {
      mln_precondition(!empty());
      auto k = h(m_ima[v]);
      while (m_hq[--k] == npos)
	;
      mln_assertion(k < h(m_ima[v]));

      auto tmp = m_hq[k];
      m_hq[k] = m_next[tmp];
      --m_size;
      return tmp;
    }


    template <typename I, typename Compare>
    inline
    typename pset<I, Compare, true>::value_type
    pset<I, Compare, true>::pop_next(const value_type& v)
    {
      mln_precondition(!empty());
      auto k = h(m_ima[v]);
      while (m_hq[k] == npos)
	++k;
      mln_assertion(not (k < h(m_ima[v])));

      auto tmp = m_hq[k];
      m_hq[k] = m_next[tmp];
      --m_size;
      return tmp;
    }



  }

}


#endif // ! MLN_MOPRHO_PSET_PSET_HPP
