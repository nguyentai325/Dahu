#ifndef LITERAL_HPP
# define LITERAL_HPP

# include <mln/core/colors.hpp>

namespace mln
{
  namespace colors
  {

    namespace literal
    {




      //struct red_t : literal_color<rgb8>;
      //using green_t = literal_color<rgb8, rgb8{0,255,0}>;
      //using blue_t = literal_color<rgb8, rgb8{0,0,255}>;
      namespace
      {
        constexpr rgb8 black{0,0,0};
        constexpr rgb8 red{255,0,0};
        constexpr rgb8 green{0,255,0};
        constexpr rgb8 blue{0,0,255};
        constexpr rgb8 white{255,255,255};
        constexpr rgb8 cyan{0,255,255};
        constexpr rgb8 yellow{255,255,0};
        constexpr rgb8 magenta{255,0,255};
      }

    }
  }
}


#endif // ! LITERAL_HPP
