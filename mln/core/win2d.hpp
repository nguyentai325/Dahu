#ifndef MLN_CORE_WIN2D_HPP
# define MLN_CORE_WIN2D_HPP

# include <mln/core/neighborhood/dyn_neighborhood.hpp>
# include <mln/core/point.hpp>
# include <mln/core/domain/box.hpp>
# include <mln/core/assert.hpp>

///
/// \file
///

namespace mln
{

  ///
  /// \brief Define a dynamic rectangular window
  ///
  struct rect2d;


  /// \defgroup Free functions
  /// \{
  rect2d make_rectangle2d(unsigned height, unsigned width);
  rect2d make_rectangle2d(unsigned height, unsigned width, point2d center);
  /// \}

  namespace {

    struct winc4_t : dyn_neighborhood_base<std::array<point2d, 5>, constant_neighborhood_tag,  winc4_t>
    {
      static const int static_size = 5;
      static constexpr std::array<point2d, 5> dpoints = {{ {-1,0}, {0,-1}, {0,0}, {0,1}, {1,0} }};
    };

    constexpr std::array<point2d, 5> winc4_t::dpoints;

    static const winc4_t winc4 {};

    struct winc8_t : dyn_neighborhood_base< std::array<point2d, 9>,  constant_neighborhood_tag, winc8_t >
    {
      static const int static_size = 9;
      static constexpr std::array<point2d, 9> dpoints = {{ {-1,-1}, {-1,0}, {-1,1},
							   {0, -1}, {0,0},  {0, 1},
							   {1,-1},  {1,1},  {1,0} }};
    };

    constexpr std::array<point2d, 9> winc8_t::dpoints;

    static const winc8_t winc8 {};

    struct winc2_v_t : dyn_neighborhood_base< std::array<point2d, 3>,  constant_neighborhood_tag, winc2_v_t >
    {
      static const int static_size = 3;
      static constexpr std::array<point2d, 3> dpoints = {{ {-1,0}, {0, 0}, {1,0} }};
    };

    constexpr std::array<point2d, 3> winc2_v_t::dpoints;

    static const winc2_v_t winc2_v {};



    struct winc2_h_t : dyn_neighborhood_base< std::array<point2d, 3>,  constant_neighborhood_tag, winc2_h_t >
    {
      static const int static_size = 3;
      static constexpr std::array<point2d, 3> dpoints = {{ {0,-1}, {0, 0}, {0,1} }};
    };

    constexpr std::array<point2d, 3> winc2_h_t::dpoints;

    static const winc2_h_t winc2_h {};

  }

  /**************************/
  /***  Implementation     **/
  /**************************/

  struct rect2d :
    dyn_neighborhood_base<box2d, dynamic_neighborhood_tag, rect2d>
  {
    typedef std::true_type                                is_incremental;
    typedef std::false_type                               is_separable;

    typedef dyn_neighborhood<box2d, dynamic_neighborhood_tag> dec_type;
    typedef dyn_neighborhood<box2d, dynamic_neighborhood_tag> inc_type;

    rect2d() = default;

    rect2d(const box2d& r)
      : dyn_neighborhood_base<box2d, dynamic_neighborhood_tag, rect2d>(r)
    {
    }

    inc_type inc() const
    {
      box2d b = this->dpoints;
      b.pmin[1] = b.pmax[1] - 1;
      return b;
    }

    dec_type dec() const
    {
      box2d b = this->dpoints;
      b.pmin[1] -= 1;
      b.pmax[1] = b.pmin[1]+1;
      return b;
    }

  };

  inline
  rect2d make_rectangle2d(unsigned height, unsigned width)
  {
    mln_precondition(height % 2 == 1);
    mln_precondition(width % 2 == 1);
    int h = height / 2;
    int w = width / 2;
    box2d b = {point2d(-h,-w), point2d(h+1,w+1)};
    return b;
  }

  inline
  rect2d make_rectangle2d(unsigned height, unsigned width, point2d center)
  {
    mln_precondition(height % 2 == 1);
    mln_precondition(width % 2 == 1);
    unsigned h = height / 2;
    unsigned w = width / 2;
    point2d uleft = center;
    point2d lright = center;
    uleft[0] -= h;
    uleft[1] -= w;
    lright[0] += h+1;
    lright[1] += w+1;
    box2d b = {uleft, lright};
    return b;
  }

}

#endif // !  MLN_CORE_WIN2D_HPP
