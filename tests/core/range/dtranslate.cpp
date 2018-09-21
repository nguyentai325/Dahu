#include <mln/core/domain/dtranslate.hpp>
#include <mln/core/domain/box.hpp>

int main()
{
  using namespace mln;

  box2d b = {{1,2}, {10,10}};
  point2d p = {3,4};


  for (auto v: dtranslate(b, p))
    std::cout << v << std::endl;
}

