#ifndef MLN_CORE_NEIGHBORHOOD_PIXEL_RANGE_HPP
# define MLN_CORE_NEIGHBORHOOD_PIXEL_RANGE_HPP


# include <mln/core/iterator/sliding_win_pixter.hpp>
# include <mln/core/std/array.hpp>

namespace mln
{

  template <typename Pixel, typename SiteSet>
  struct neighborhood_pixel_range;


  /******************************************/
  /****          Implementation          ****/
  /******************************************/


  template <typename Pixel, typename SiteSet>
  struct neighborhood_pixel_range
  {
    typedef sliding_win_pixter<SiteSet, Pixel> iterator;
    typedef iterator const_iterator;

    neighborhood_pixel_range(const SiteSet& s, const Pixel& pix)
      : s_ (s), pix_ (pix)
    {
    }

    iterator iter() const
    {
      return iterator(s_, pix_);
    }

  private:
    const SiteSet& s_;
    const Pixel& pix_;
  };


} // end of namespace mln


#endif //!MLN_CORE_NEIGHBORHOOD_PIXEL_RANGE_HPP
