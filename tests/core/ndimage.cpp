#include <mln/core/ndimage.hpp>
#include <mln/core/grays.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm_ext.hpp>

int main()
{
  using namespace mln;
  image2d<uint8> ima( (box2d) {{2,5}, {4,9}} );

  boost::iota(ima.values(), 1);

  std::cout << "== Iteration on values == " << std::endl;

  std::cout << "FWD: ";
  for (int v: ima.values())
    std::cout << v << ",";
  std::cout << std::endl;

  std::cout << "BWD: ";
  for (int v: boost::adaptors::reverse(ima.values()))
    std::cout << v << ",";
  std::cout << std::endl;


  std::cout << "== Iteration on pixels == " << std::endl;

  std::cout << "FWD: ";
  for (auto v: ima.pixels())
    std::cout << v.point() << ":" << (int)v.val() << ",";
  std::cout << std::endl;

  std::cout << "BWD: ";
  for (auto v: boost::adaptors::reverse(ima.pixels()))
    std::cout << v.point() << ":" << (int)v.val() << ",";
  std::cout << std::endl;


}
