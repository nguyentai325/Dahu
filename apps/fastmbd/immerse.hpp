#ifndef IMMERSE_HPP
#define IMMERSE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/core/image/sub_image.hpp>
#include "range.hpp"



namespace mln
{

    bool is_border (point2d p, const image2d<uint8_t>& ima)
    {
        unsigned height = ima.nrows();
        unsigned width = ima.ncols();

        if ((p[0] >= 0 and p[0]  <= 3)   or (p[0] <= height -1 and p[0] >= height - 4) or (p[1] >= 0 and p[1] <= 3)   or (p[1]  <= width - 1 and p[1] >= width - 4))
            return true;
        else
            return false;
    }

    image2d< range<uint8_t> >
    add_kfaces(const image2d<uint8_t>& ima)
    {

      box2d b = ima.domain();
      b.pmin = b.pmin * 2;
      b.pmax = b.pmax * 2 - 1;

      //   0 1
      // 0 a b
      // 1 c d

      // ->

      //  -1 0 1 2 3
      //-1 + - + - +
      // 0 | a | b |
      // 1 + - + - +
      // 2 | c | d |
      // 3 + - + - +

      typedef range<uint8_t> R;
      typedef point2d      P;
      image2d<R> out(b);

      mln_foreach(point2d p, ima.domain())
        {

          uint8_t a = ima.at(p),
            b = ima.at(p + P{0,1}),
            c = ima.at(p + P{1,0}),
            d = ima.at(p + P{1,1});

          uint8_t min1 = std::min(a,b), min2 = std::min(a,c);
          uint8_t max1 = std::max(a,b), max2 = std::max(a,c);
          uint8_t min3 = std::min(d, std::min(c, min1));
          uint8_t max3 = std::max(d, std::max(c, max1));

          point2d q = 2 * p;
          out.at(q) = ima.at(p);

          if (is_border(q + P{0,1}, ima))
              out.at(q + P{0,1}) = min1;
          else
              out.at(q + P{0,1}) = R{min1, max1};

          if (is_border(q + P{1,0}, ima))
              out.at(q + P{1,0}) = min2;
          else
              out.at(q + P{1,0}) = R{min2, max2};

          out.at(q + P{1,1}) = R{min3, max3};

        }

      return out;
    }


    image2d< range<uint8_t> >
    immerse(const image2d<uint8_t>& input)
    {
      // image2d<value::int_u8> temp = add_border(input);
      // temp = subdivide(temp);
      // return add_kfaces(temp);
      return add_kfaces(input);
    }
}

#endif // IMMERSE_HPP
