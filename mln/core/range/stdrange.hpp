#ifndef STDRANGE_HPP
# define STDRANGE_HPP

namespace mln
{

  // Mixin pattern
  // Provide mln style iterator for a std range
  //
  template <typename StandardRange>
  struct stdrange : StandardRange
  {
  private:
    typedef StandardRange base;

  public:
    stdrange(const StandardRange& r) :
    base(r)
    {
    }

    stdrange(StandardRange&& r) :
    base (std::move(r))
    {
    }

    stdrange(const stdrange& r) :
    base (static_cast<const base&>(r))
    {
    }

    stdrange(stdrange&& r) :
    base (static_cast<base&&>(r))
    {
    }

    typedef stditerator<typename StandardRange::iterator> iterator;
    typedef stditerator<typename StandardRange::const_iterator> const_iterator;

    iterator iter()
    {
      iterator (base::begin(), base::end());
    }

    const_iterator iter() const
    {
      const_iterator (base::begin(), base::end());
    }
  };

  // Free functions
  template <typename StandardRange>
  stdrange<StandardRange>
  make_stdrange(StandardRange&& x)
  {
    return stdrange<StandardRange>(std::forward<StandardRange>(x));
  }

  template <typename StandardRange>
  stdrange<StandardRange>&
  make_stdrange(StandardRange& x)
  {
    return static_cast<stdrange<StandardRange>&>(x);
  }
}


#endif // ! STDRANGE_HPP
