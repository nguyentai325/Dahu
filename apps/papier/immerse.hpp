#ifndef IMMERSE_HPP
#define IMMERSE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/core/image/sub_image.hpp>
#include "range.hpp"



namespace mln
{
    typedef std::array<range<uint8>, 3> range_vec;

    template <typename T>
    bool is_border (point2d p, const image2d<T>& ima)
    {
        unsigned height = ima.nrows();
        unsigned width = ima.ncols();

        if ((p[0] == 0)   or (p[0] == height - 1) or (p[1] == 0)   or (p[1]  == width - 1))
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

          if (is_border(q + P{1,1}, ima))
              out.at(q + P{1,1}) = min3;
          else
              out.at(q + P{1,1}) = R{min3, max3};

        }

      return out;
    }


    image2d< range_vec >
    add_kfaces_color(const image2d<rgb8>& ima)
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
      image2d<range_vec> out(b);

      mln_foreach(point2d p, ima.domain())
        {

          rgb8 a = ima.at(p),
            b = ima.at(p + P{0,1}),
            c = ima.at(p + P{1,0}),
            d = ima.at(p + P{1,1});
          for (int i = 0 ; i < 3; i++)
          {
              uint8_t min1 = std::min(a[i],b[i]), min2 = std::min(a[i],c[i]);
              uint8_t max1 = std::max(a[i],b[i]), max2 = std::max(a[i],c[i]);
              uint8_t min3 = std::min(d[i], std::min(c[i], min1));
              uint8_t max3 = std::max(d[i], std::max(c[i], max1));

              point2d q = 2 * p;
              out.at(q)[i] = ima.at(p)[i];

              if (is_border(q + P{0,1}, ima))
                  out.at(q + P{0,1})[i] = min1;
              else
                  out.at(q + P{0,1})[i] = R{min1, max1};

              if (is_border(q + P{1,0}, ima))
                  out.at(q + P{1,0})[i] = min2;
              else
                  out.at(q + P{1,0})[i] = R{min2, max2};

              if (is_border(q + P{1,1}, ima))
                  out.at(q + P{1,1})[i] = min3;
              else
                  out.at(q + P{1,1})[i] = R{min3, max3};
          }

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

    image2d< range_vec >
    immerse_color(const image2d<rgb8>& input)
    {
      return add_kfaces_color(input);
    }
}

#endif // IMMERSE_HPP
