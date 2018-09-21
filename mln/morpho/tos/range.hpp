#ifndef RANGE_HH
# define RANGE_HH

# include <cstdlib>
# include <iostream>


namespace mln
{

  namespace tos
  {


    // range


    template <typename T> 
    struct range
    {
      typedef T value;

      T lower, upper;

      range()
      {
      }

      template <typename T_>
      range(T_ v)
      {
	lower = upper = v;
      }

      range(T l, T u)  // permissive to args ordering
      {
	if (l <= u)
	  {
	    this->lower = l;
	    this->upper = u;
	  }
	else
	  {
	    this->lower = u;
	    this->upper = l;
	  }
      }

      range& operator=(const range& rhs)
      {
	lower = rhs.lower;
	upper = rhs.upper;
	return *this;
      }

      template <typename T_>
      range& operator=(T_ v)
      {
	lower = upper = v;
	return *this;
      }

      bool has(const T& v) const
      {
	mln_precondition(upper >= lower);
	return lower <= v && v <= upper;
      }

      void
      span_to(T v)
      {
	mln_precondition(upper >= lower);
	if (this->has(v))
	  return;
	if (v < lower)
	  lower = v;
	if (v > upper)
	  upper = v;
      }

      unsigned index_of(const T& v) const
      {
	mln_precondition(upper >= lower);
	if (! this->has(v))
	  std::abort();
	if (v < lower or v > upper)
	  std::abort();
	return v - lower;
      }

      T value_at(unsigned i) const
      {
	mln_precondition(upper >= lower);
	T v;
	set_enc(v, lower + i);
	if (! this->has(v))
	  std::abort();
	return v;
      }

      unsigned nelements() const
      {
	mln_precondition(upper >= lower);
	return upper - lower + 1;
      }

      T length() const
      {
	mln_precondition(upper >= lower);
	return upper - lower;
      }

      bool is_degenerated() const
      {
	mln_precondition(upper >= lower);
	return upper == lower;
      }

      T single_value() const
      {
	mln_precondition(upper >= lower);
	if (! is_degenerated())
      	  std::abort();
      	return upper;
      }

      operator T() const
      {
	mln_precondition(upper >= lower);
      	if (! is_degenerated())
      	  std::abort();
      	return upper;
      }

    };



    //  comparison

    template <typename T> 
    bool
    operator==(const range<T>& lhs, const range<T>& rhs)
    {
      return lhs.lower == rhs.lower && lhs.upper == rhs.upper;
    }

    template <typename T> 
    bool
    operator!=(const range<T>& lhs, const range<T>& rhs)
    {
      return ! (lhs == rhs);
    }


    // deactivation of ordering related operators

    template <typename T> 
    void operator<(const range<T>&, const range<T>&);

    template <typename T> 
    void operator<=(const range<T>&, const range<T>&);

    template <typename T> 
    void operator>(const range<T>&, const range<T>&);

    template <typename T> 
    void operator>=(const range<T>&, const range<T>&);



    //  set ops

    template <typename T>
    bool
    are_adjacent(const range<T>& r1, const range<T>& r2)
    {
      return span_op(r1, r2).nelements() == r1.nelements() + r2.nelements();
    }

    template <typename T>
    bool
    do_intersect(const range<T>& r1, const range<T>& r2)
    {
      return span_op(r1, r2).nelements() < r1.nelements() + r2.nelements();
    }

    template <typename T>
    range<T>
    inter(const range<T>& r1, const range<T>& r2)
    {
      if (! do_intersect(r1, r2))
	std::abort();
      return range<T>(std::max(r1.lower, r2.lower),
		      std::min(r1.upper, r2.upper));
    }




    //  span operator

    template <typename T> 
    range<T>
    span_op(const range<T>& r1, const range<T>& r2)
    {
      return range<T>(std::min(r1.lower, r2.lower),
		      std::max(r1.upper, r2.upper));
    }


    //  op<<

    template <typename T> 
    std::ostream&
    operator<<(std::ostream& ostr, const range<T>& i)
    {
      if (i.is_degenerated())
	return ostr << i.lower;
      else
	return ostr << '[' << i.lower << ',' << i.upper << ']';
    }




  } // end of namespace mln::tos

} // end of namespace mln


#endif // ndef RANGE_HH
