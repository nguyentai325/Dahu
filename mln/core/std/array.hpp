#ifndef ARRAY_HPP
# define ARRAY_HPP

# include <array>
# include <mln/core/iterator/stditerator.hpp>

namespace mln
{

  template <typename T, std::size_t N>
  struct array
  {
    typedef T				value_type;
    typedef std::size_t			size_type;
    typedef std::ptrdiff_t		distance_type;
    typedef T&				reference;
    typedef const T&			const_reference;
    typedef T*				pointer;
    typedef const T*			const_pointer;
    typedef stditerator<T*>		iterator;
    typedef stditerator<const T*>	const_iterator;
    typedef stditerator< std::reverse_iterator<T*> >	   reverse_iterator;
    typedef stditerator< std::reverse_iterator<const T*> > const_reverse_iterator;

    reference at(size_type pos)
    {
      return x_.at(pos);
    }

    const_reference at(size_type pos) const
    {
      return x_.at(pos);
    }

    reference operator[](size_type pos) noexcept
    {
      return x_[pos];
    }

    const_reference operator[](size_type pos) const noexcept
    {
      return x_[pos];
    }

    reference front() noexcept
    {
      return x_.front();
    }

    const_reference front() const noexcept
    {
      return x_.front();
    }

    reference back() noexcept
    {
      return x_.back();
    }

    const_reference back() const noexcept
    {
      return x_.back();
    }

    pointer data() noexcept
    {
      return x_.data();
    }

    const_pointer data() const noexcept
    {
      return x_.data();
    }

    constexpr bool empty() noexcept
    {
      return x_.empty();
    }

    constexpr size_type size() noexcept
    {
      return x_.size();
    }

    constexpr size_type max_size() noexcept
    {
      return x_.max_size();
    }

    size_type fill(const T& value)
    {
      return x_.fill(value);
    }

    void swap(const T& value) noexcept(noexcept(swap(std::declval<T&>(), std::declval<T&>())))
    {
      return x_.swap(value);
    }


    operator std::array<T, N>& () noexcept
    {
      return x_;
    }

    operator const std::array<T, N>& () const noexcept
    {
      return x_;
    }

    iterator iter() noexcept
    {
      return iterator(x_.begin(), x_.end());
    }

    const_iterator iter() const noexcept
    {
      return const_iterator(x_.begin(), x_.end());
    }

    iterator riter() noexcept
    {
      return iterator(x_.rbegin(), x_.rend());
    }

    const_iterator riter() const noexcept
    {
      return const_iterator(x_.rbegin(), x_.rend());
    }

    std::array<T, N> x_;
  };

}

#endif // ! ARRAY_HPP
