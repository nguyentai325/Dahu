#ifndef MLN_CORE_ITERATOR_SLIDING_WIN_PITER_HPP
# define MLN_CORE_ITERATOR_SLIDING_WIN_PITER_HPP

# include <mln/core/iterator/iterator_base.hpp>

namespace mln
{

  template <typename SiteSet>
  struct sliding_win_piter;


  /// Specialization for static sized domain.
  template <typename Point, size_t N>
  struct sliding_win_piter< mln::array<Point, N> >;


  /******************************************/
  /****          Implementation          ****/
  /******************************************/



  /******************************************/
  /****          Specialization          ****/
  /******************************************/


  template <typename Point, size_t N>
  struct sliding_win_piter< mln::array<Point, N> > :
    iterator_base< sliding_win_piter< mln::array<Point, N> >,
                   const Point >
  {
    sliding_win_piter() = default;
    sliding_win_piter(const mln::array<Point, N>& siteset, const Point& p)
      : arr_(siteset), p_ (&p)
    {
    }

      void init() { i_ = 0; cur_ = *p_ + arr_[0]; }
    void next() { cur_ = *p_ + arr_[++i_]; }
    bool finished() const { return i_ == N; }
    const Point& dereference() const { return cur_; }

    // Specific to sliding iterator
    void center(const Point& p) { p_ = &p; }

  private:
    mln::array<Point, N> arr_;
    const Point* p_;
    int i_;
    Point cur_;
  };



} // end of namespace mln

#endif //!MLN_CORE_ITERATOR_SLIDING_WIN_PITER_HPP
