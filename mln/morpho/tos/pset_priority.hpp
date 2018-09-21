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

#ifndef MLN_MOPRHO_TOS_PSET_PRIORITY_HPP
# define MLN_MOPRHO_TOS_PSET_PRIORITY_HPP

# include <set>

namespace mln
{

  namespace morpho
  {

    //bool Enable = has_indexer<typename I::value_type, Compare>::value  >
    template <typename I,
	      typename Compare = std::less<typename I::value_type>,
	      bool Enable = false>
    struct pset_priority;


    template <typename I, typename Compare>
    struct pset_priority<I, Compare, false>
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

      pset_priority(const I& ima, const Compare& cmp = Compare() );

      void insert(const value_type& v);
      bool empty() const;


      /// \brief returns the next point \p k such that level(p) <= level(k)
      /// If there is point \p k such that level(p) == level(k) then k is returned
      /// Otherwise it returns the point at highest priority (rend())
      bool		has_next(const value_type& v) const;
      value_type	pop_next(const value_type& p);

      /// \brief returns the next point \p k such that level(p) < level(k)
      bool		has_previous(const value_type& v) const;
      value_type	pop_previous(const value_type& p);

    private:
      std::multiset<value_type, cmp_t>		m_set;
    };

    /*************************/
    /**  Implementation     **/
    /*************************/



    template <typename I, typename Compare>
    pset_priority<I, Compare, false>::pset_priority(const I& ima, const Compare& cmp)
      : m_set(cmp_t {ima, cmp} )
    {
    }

    template <typename I, typename Compare>
    void
    pset_priority<I, Compare, false>::insert(const value_type& v)
    {
      m_set.insert(v);
    }

    template <typename I, typename Compare>
    bool
    pset_priority<I, Compare, false>::empty() const
    {
      return m_set.empty();
    }

    template <typename I, typename Compare>
    bool
    pset_priority<I, Compare, false>::has_previous(const value_type& p) const
    {
      mln_precondition(!empty());
      return m_set.key_comp()(*(m_set.begin()), p);
    }

    template <typename I, typename Compare>
    bool
    pset_priority<I, Compare, false>::has_next(const value_type& p) const
    {
      mln_precondition(!empty());
      return !m_set.key_comp()(*(m_set.rbegin()), p);
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset_priority<I, Compare, false>::pop_previous(const value_type& p)
    {
      mln_precondition(has_previous(p));
      auto x = m_set.lower_bound(p);
      --x;
      mln_assertion(x != m_set.end());
      value_type v = *x;
      m_set.erase(x);
      return v;
    }

    template <typename I, typename Compare>
    typename I::size_type
    pset_priority<I, Compare, false>::pop_next(const value_type& p)
    {
      mln_precondition(has_next(p));
      auto x = m_set.lower_bound(p); // level(p) <= level(x)
      mln_assertion(x != m_set.end());

      if (not m_set.key_comp()(p, *x)) // then level(x) == level(p)
	;
      else // level(p) < level(x) => get the highest priority 
	x = --m_set.end();

      value_type v = *x;
      m_set.erase(x);
      return v;
    }

  }

}


#endif // ! MLN_MOPRHO_TOS_PSET_PRIORITY_HPP
