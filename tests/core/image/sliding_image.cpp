#include <mln/core/image/sliding_image.hpp>
#include <mln/core/ndimage.hpp>
#include <mln/core/grays.hpp>
#include <mln/io/imprint.hpp>
#include <array>

int main()
{
  using namespace mln;
  typedef std::array<point2d, 3> S;
  S delta = {{ {{-1,-1}}, {{0,0}}, {{1,1}} }};

  image2d<uint8> ima(5, 5);
  sliding_image<image2d<uint8>, S> f(ima, delta);

  f.center(point2d{{1,1}});
  io::imprint(f);

  f.center(point2d{{3,2}});
  io::imprint(f);
}
