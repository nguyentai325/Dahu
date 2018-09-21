
#include <vector>
#include <iostream>
#include <mln/core/domain/box.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/timer.hpp>

int main()
{
  std::vector<int> s = {1,2,3,4,5};


  std::cout << "Direct: ";
  for(auto v: s)
    std::cout << v << ", ";
  std::cout << std::endl;


  std::cout << "Reverse: ";
  for(auto v: boost::adaptors::reverse(s) )
    std::cout << v << ", ";
  std::cout << std::endl;

  mln::box2d b = {{2,3}, {4,7}};
    std::cout << "Direct: ";
  for(auto p: b)
    std::cout << p << ", ";
  std::cout << std::endl;


  std::cout << "Reverse: " << std::endl;
  for(auto p: boost::adaptors::reverse(b) )
    std::cout << p << ", ";
  std::cout << std::endl;

}

