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

#ifndef INDEX_HPP
# define INDEX_HPP

# include <mln/core/value/value_traits.hpp>
//# include <boost/operators.hpp>
# include <functional>

namespace mln
{

  template <typename T, typename StrictWeakOrdering>
  struct Index;


  template <typename T>
  struct Index<T, std::less<T> >
  {
    Index() = default;
    explicit Index(T v) : x (v) {}

    // operator long () const { return x; }
    // operator unsigned long () const { return x; }
    operator int () const { return x; }
    // operator unsigned () const { return x; }
    // operator T () const { return x; }

    Index& operator++ ()    { ++x; return *this; }
    Index& operator-- ()    { --x; return *this; }
    Index  operator++ (int) { return Index(x++); }
    Index  operator-- (int) { return Index(x--); }
    bool operator== (const Index& other) const { return x == other.x; }
    bool operator!= (const Index& other) const { return x != other.x; }
    bool operator< (const Index& other)  const { return x < other.x; }
    bool operator> (const Index& other)  const { return x > other.x; }
    bool operator<= (const Index& other) const { return x <= other.x; }
    bool operator>= (const Index& other) const { return x >= other.x; }
  private:
    T x;
  };

  template <typename T>
  struct Index<T, std::greater<T> >
  {
    Index() = default;
    explicit Index(T v) : x (v) {}

    // operator long () const { return x; }
    // operator unsigned long () const { return x; }
    operator int () const { return x; }
    // operator unsigned () const { return x; }
    //operator T () const { return x; }


    Index& operator++ ()    { --x; return *this; }
    Index& operator-- ()    { ++x; return *this; }
    Index  operator++ (int) { return Index(x--); }
    Index  operator-- (int) { return Index(x++); }
    bool operator== (const Index& other) const { return x == other.x; }
    bool operator!= (const Index& other) const { return x != other.x; }
    bool operator< (const Index& other)  const { return x > other.x; }
    bool operator> (const Index& other)  const { return x < other.x; }
    bool operator<= (const Index& other) const { return x >= other.x; }
    bool operator>= (const Index& other) const { return x <= other.x; }
  private:
    T x;
  };



  template <typename T>
  struct value_traits< Index<T, std::less<T> >, std::less< Index<T, std::less<T> > >, void>
  {
  private:
    typedef Index<T, std::less<T> > index_t;
  public:
    static constexpr unsigned quant = value_traits<T>::quant;
    static constexpr index_t min() { return index_t(0); }
    static constexpr index_t max() { return index_t(value_traits<T>::max()); }
    static constexpr index_t inf() { return min(); }
    static constexpr index_t sup() { return max(); }
  };

  template <typename T>
  struct value_traits< Index<T, std::greater<T> >, std::less< Index<T, std::greater<T> > >, void>
  {
  private:
    typedef Index<T, std::greater<T> > index_t;
  public:
    static constexpr unsigned quant = value_traits<T>::quant;
    static constexpr index_t min() { return index_t(value_traits<T>::max()); }
    static constexpr index_t max() { return index_t(0); }
    static constexpr index_t inf() { return min(); }
    static constexpr index_t sup() { return max(); }
  };

  template <typename T>
  struct value_traits< Index<T, std::less<T> >, productorder_less<Index<T, std::less<T> > >, void>
    : value_traits< Index<T, std::less<T> >, std::less< Index<T, std::less<T> > >, void>
  {
  };

  template <typename T>
  struct value_traits< Index<T, std::greater<T> >, productorder_less<Index<T, std::greater<T> > >, void>
    : value_traits< Index<T, std::greater<T> >, std::less< Index<T, std::greater<T> > >, void>
  {
  };

}


#endif // ! INDEX_HPP

