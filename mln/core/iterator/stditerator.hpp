#ifndef STDITERATOR_HPP
# define STDITERATOR_HPP

# include <mln/core/iterator/iterator_base.hpp>
# include <type_traits>
# include <iterator>

namespace mln
{

  template <typename Iterator>
  struct stditerator : iterator_base< stditerator<Iterator>,
				      typename std::iterator_traits<Iterator>::value_type,
				      typename std::iterator_traits<Iterator>::reference >
  {
    typedef typename std::iterator_traits<Iterator>::reference reference;

    stditerator()
    {
    }

    stditerator(const Iterator& begin, const Iterator& end)
      : begin_(begin), end_( end)
    {
    }

    template <typename Iterator2>
    stditerator(const stditerator<Iterator2>& other,
		typename std::enable_if<std::is_convertible<Iterator2, Iterator>::value>::type* = NULL)
      : cur_(other.cur_), begin_ (other.begin_), end_(other.end_)
    {
    }

    void init()
    {
      cur_ = begin_;
    }

    void next()
    {
      ++cur_;
    }

    bool finished() const
    {
      return cur_ == end_;
    }

    reference dereference() const
    {
      return *cur_;
    }

  private:
    template <typename> friend struct stditerator;

    Iterator cur_;
    Iterator begin_;
    Iterator end_;
  };

}

#endif // ! STDITERATOR_HPP
