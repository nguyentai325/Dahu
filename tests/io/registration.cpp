#include <mln/io/registration.hpp>
#include <iostream>

int main()
{
  {
    std::map<int, const mln::io::type_info_base*>::const_iterator x =
      mln::io::registration_map_.begin();

    for (x; x !=  mln::io::registration_map_.end(); ++x)
      std::cout << x->first << ":" << x->second->name() << std::endl;
  }

  mln::io::type_info<mln::uint8> t;
  std::cout << t.id() << ":" << t.name() << std::endl;

  {
    std::map<int, const mln::io::type_info_base*>::const_iterator x =
      mln::io::registration_map_.begin();

    for (x; x !=  mln::io::registration_map_.end(); ++x)
      std::cout << x->first << ":" << x->second->name() << std::endl;
  }
}
