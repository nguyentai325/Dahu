#include <mln/core/ndimage.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/macros.hpp>
#include <mln/io/imprint.hpp>

#include <boost/range/algorithm_ext/iota.hpp>

int main()
{
  using namespace mln;

  image2d<uint8> ima(5, 5);
  boost::iota(ima.values(), 1);


  io::imprint(ima);

  // Test point iterator
  std::cout << "Testing value iterator" << std::endl;
  forallp(p, ima) {
    std::cout << p << ":" << (int)ima(p) << std::endl;
  }


  // Test pixel iterator
  std::cout << "Testing value iterator" << std::endl;
  forallv(v, ima) {
    std::cout << (int)v << std::endl;
  }


  // Test pixel iterator
  std::cout << "Testing pixel iterator" << std::endl;
  forall(x, ima) {
    std::cout << x.point() << " : " << (int)x.val() << std::endl;
  }


}
