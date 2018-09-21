#ifndef MLN_CORE_NEIGHB2D_HPP
# define MLN_CORE_NEIGHB2D_HPP

//# include <mln/core/neighborhood.hpp>
//# include <mln/core/std/array.hpp>
# include <array>

# include <mln/core/neighborhood/dyn_neighborhood.hpp>
# include <mln/core/point.hpp>


namespace mln {

  namespace {

    struct c4_t : dyn_neighborhood_base<std::array<point2d, 4>, constant_neighborhood_tag,  c4_t>
    {
      static const int static_size = 4;
      static constexpr std::array<point2d, 4> dpoints = {{ {-1,0}, {0,-1}, {0,1}, {1,0} }};
    };


    struct c8_t : dyn_neighborhood_base< std::array<point2d, 8>,  constant_neighborhood_tag, c8_t >
    {
      static const int static_size = 8;
      static constexpr std::array<point2d, 8> dpoints = {{ {-1,-1}, {-1,0}, {-1,1},
							   {0, -1},         {0, 1},
							   {1,-1},  {1,0},  {1,1} }};
    };


    struct c2_v_t : dyn_neighborhood_base< std::array<point2d, 2>,  constant_neighborhood_tag, c2_v_t >
    {
      static const int static_size = 2;
      static constexpr std::array<point2d, 2> dpoints = {{ {-1,0}, {1,0} }};
    };


    struct c2_h_t : dyn_neighborhood_base< std::array<point2d, 2>,  constant_neighborhood_tag, c2_h_t >
    {
      static const int static_size = 2;
      static constexpr std::array<point2d, 2> dpoints = {{ {0,-1}, {0,1} }};
    };


  constexpr std::array<point2d, 4> c4_t::dpoints;
  constexpr std::array<point2d, 8> c8_t::dpoints;
  constexpr std::array<point2d, 2> c2_v_t::dpoints;
  constexpr std::array<point2d, 2> c2_h_t::dpoints;
}


  static const c4_t c4 {};
  static const c8_t c8 {};
  static const c2_v_t c2_v {};
  static const c2_h_t c2_h {};

} // end of namespace mln

#endif //!MLN_CORE_NEIGHB2D_HPP
