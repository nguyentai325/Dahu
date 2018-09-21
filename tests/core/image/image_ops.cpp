#include <mln/core/image/image2d.hpp>
#include <mln/core/image/image_ops.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/algorithm/fill.hpp>

#include <mln/io/imprint.hpp>
#include <boost/tuple/tuple_io.hpp>
#define BOOST_TEST_MODULE Core
#include <boost/test/unit_test.hpp>


mln::image2d<int>
make_image()
{
  mln::image2d<int> x(5, 5);
  mln::iota(x, 0);
  return x;
}

struct rgb
{
  int r,g,b;

  bool
  operator== (const rgb& other) const
  {
    return r == other.r and g == other.g and b == other.b;
  }
};

struct red : std::unary_function<rgb&, int&>
{

  int& operator() (rgb& x) const
  {
    return x.r;
  }

  const int&
  operator() (const rgb& x) const
  {
    return x.r;
  }
};

std::ostream&
operator<< (std::ostream& ss, const rgb& x)
{
  return ss << boost::make_tuple(x.r, x.g, x.b);
}


BOOST_AUTO_TEST_CASE(LValueOperator)
{
  using namespace mln;
  image2d<rgb> ima(5,5);

  rgb zero = {0,0,0};
  rgb douze = {12,0,0};
  fill(ima, zero);

  auto x = make_unary_image_expr(ima, red ());
  fill(x, 12);

  BOOST_CHECK( all(ima == douze) );
}

BOOST_AUTO_TEST_CASE(Operators)
 {
   using namespace mln;


  image2d<int> ima(5,5);
  image2d<int> ref(5,5);

  iota(ima, 0);
  int i = 0;

  mln_viter(v, ref);
  mln_forall(v)
    *v = i--;

  BOOST_CHECK( all(-ima == ref) );
 }



BOOST_AUTO_TEST_CASE(MixedOperator)
{
  using namespace mln;

  image2d<char>	  x(5,5);
  image2d<short>  y(5,5);

  iota(x, 0);
  iota(y, 0);

  BOOST_CHECK( (std::is_same<typename decltype(x + x)::value_type, char> ()) );
  BOOST_CHECK( (std::is_same<typename decltype(x + y)::value_type, typename std::common_type<char, short>::type> ()) );


  BOOST_CHECK( all((x + y) == (2*y)) );
}


//   // Const: Lvalue + Lvalue
//   io::imprint(ima + ima);

//   auto x = make_image() + make_image();
//   io::imprint(x * x / 2);


//   // Const: Rvalue + Lvalue
//   io::imprint((ima + ima) - (ima));
// }

