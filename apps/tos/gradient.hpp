#ifndef GRADIENT_HPP
# define GRADIENT_HPP

# include <mln/core/canvas/accfpfn.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/win2d.hpp>


namespace mln
{

  template <typename T>
  struct Minmax : Accumulator< Minmax<T> >
  {
    typedef T argument_type;
    typedef std::pair<T,T> result_type;

    void init() {
      min = value_traits<T>::max();
      max = value_traits<T>::min();
    }

    void take(const T& v) {
      if (v < min)
	min = v;
      if (v > max)
	max = v;
    }


    T min, max;
  };




  template <typename T>
  image2d<T>
  gradient(const image2d<T>& ima, int size)
  {
    mln_precondition( size > 1 );

    rect2d win = make_rectangle2d(size, size);
    return accfpfn(ima, win, true, Minmax<T> (), [] (const Minmax<T>& a) -> T {
	return a.max - a.min; });
  }


}

#endif // ! GRADIENT_HPP
