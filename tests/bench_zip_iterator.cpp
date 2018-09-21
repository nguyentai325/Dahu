#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/iota.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/core/algorithm/equal.hpp>
#include <mln/core/extension/fill.hpp>

#include <mln/core/range/algorithm/accumulate.hpp>
#include <mln/core/range/algorithm/inner_product.hpp>

#include <mln/io/imprint.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/wrt_offset.hpp>

//#include <mln/core/forall.hpp>
//#include <mln/core/pix_range_iterator_proxy.hpp>

#include <mln/kernel/kernel.hpp>
#include <mln/accu/accumulators/min.hpp>
#include <mln/accu/accumulators/max.hpp>

#include <boost/range/numeric.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace mln;

const int& access_point(const image2d<int>& a, point2d p)
{
  const int& x = a(p);
  asm ("");
  return x;
}

const int& access_index(const image2d<int>& a, std::size_t i)
{
  const int& x = a[i];
  asm ("");
  return x;
}


double test_native(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_native(const image2d<int>& a, const image2d<int>& b)
{
  const size_t* strides = a.strides();
  const char* ptra = (const char*) &(a(point2d {0,0}));
  const char* ptrb = (const char*) &(b(point2d {0,0}));

  double v = 0;
  unsigned nrows = a.nrows(), ncols = a.ncols();
  for (unsigned i = 0; i < nrows; ++i, ptra += strides[0], ptrb += strides[0])
    {
      const int* pa = (const int*) ptra, *pb = (const int*) ptrb;
      for (unsigned j = 0; j < ncols; ++j, ++pa, ++pb)
	v += *pa * *pb;
    }

  return v;
}

double test_iterator(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_iterator(const image2d<int>& a, const image2d<int>& b)
{
  double v = range::inner_product(a.values(), b.values(), 0.0);

  return v;
}

double test_zip(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_zip(const image2d<int>& a, const image2d<int>& b)
{
  auto ima = imzip(a, b);

  double v = range::accumulate(ima.values(), 0.0, [](double v, const std::tuple<const int&, const int&>& x)
			       { return v + std::get<0>(x) * std::get<1>(x); });

  return v;
}

double test_zip_pix(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_zip_pix(const image2d<int>& a, const image2d<int>& b)
{
  auto ima = imzip(a, b);
  typedef typename decltype(ima)::const_pixel_type pixel_t;
  double v = range::accumulate(ima.pixels(), 0.0, [](double v, const pixel_t& x)
                               { return v + std::get<0>(x.val()) * std::get<1>(x.val()); });
  return v;
}


double test_for(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_for(const image2d<int>& a, const image2d<int>& b)
{
  double v = 0;
  mln_viter(x, y, a, b);

  mln_forall(x, y)
    v += *x * *y;

  return v;
}

double test_for_pixel(const image2d<int>& a, const image2d<int>& b) __attribute__ ((noinline));

double test_for_pixel(const image2d<int>& a, const image2d<int>& b)
{
  double v = 0;

  mln_pixter(x, y, a, b);
  mln_forall (x, y)
    v += x->val() * y->val();

  return v;
}

void test_dilation_native(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));

void test_dilation_native(const image2d<int>& a, image2d<int>& b)
{
  extension::fill(b, std::numeric_limits<int>::max());

  const size_t* strides = a.strides();
  const char* ptra = (const char*) &(a(point2d {0,0}));
  char* ptrb = (char*) &(b(point2d {0,0}));

  unsigned nrows = a.nrows(), ncols = a.ncols();
  auto offsets = wrt_offset(a, c8_t::dpoints);

  for (unsigned i = 0; i < nrows; ++i, ptra += strides[0], ptrb += strides[0])
    {
      const int* pa = (const int*) ptra;
      int *pb = (int*) ptrb;
      for (unsigned j = 0; j < ncols; ++j, ++pa, ++pb)
	{
	  int tmp = value_traits<int>::max();
	  for (int k = 0; k < 8; ++k)
	    tmp = std::min(tmp, *(const int*) ((const char*)pa + offsets[k]));
	  *pb = tmp;
	}
    }
}


void test_dilation_pixel(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));

void test_dilation_pixel(const image2d<int>& a, image2d<int>& b)
{
  extension::fill(b, std::numeric_limits<int>::max());

  mln_pixter(pin, pout, a, b);
  mln_iter(x, c8(pin));


  mln_forall(pin, pout)
    {
      int tmp = value_traits<int>::max();
      mln_forall(x)
	tmp = std::min(tmp, x->val());
      pout->val() = tmp;
    }
}

void test_dilation_point(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));
void test_dilation_point(const image2d<int>& a, image2d<int>& b)
{
  extension::fill(b, std::numeric_limits<int>::max());

  mln_iter(p, a.domain());
  mln_iter(n, c8(*p));

  mln_forall(p)
  {
    int tmp = value_traits<int>::max();
    mln_forall(n)
      tmp = std::min(tmp, a.at(*n));
    b(*p) = tmp;
  }
}

void test_dilation_pixel_2(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));
void test_dilation_pixel_2(const image2d<int>& a, image2d<int>& b)
{
  extension::fill(b, std::numeric_limits<int>::max());

  mln_pixter(pin, a);
  mln_viter(pout, b);

  //mln_pixter(tmp, a);
  //mln_iter(x, c8(tmp));
  mln_iter(x, c8(pin));

  //tmp.init();
  for(pin.init(), pout.init(); !pin.finished(); pin.next(), pout.next())
    {
      int tmp = value_traits<int>::max();
      mln_forall(x)
	tmp = std::min(tmp, x->val());
      *pout = tmp;
    }
}

void test_dilation_index(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));
void test_dilation_index(const image2d<int>& a, image2d<int>& b) 
{
  extension::fill(b, std::numeric_limits<int>::max());

  std::size_t idx = a.index_of_point(a.domain().pmin);
  auto w = wrt_delta_index(a, c8_t::dpoints);

  unsigned nrows = a.nrows();
  unsigned ncols = a.ncols();
  for (unsigned i = 0; i < nrows; ++i, idx += a.index_strides()[0])
    {
      for (unsigned j = 0; j < ncols; ++j)
        {
          unsigned p = idx + j;
	  int tmp = value_traits<int>::max();
          for (unsigned k = 0; k < 8; k++) {
            tmp = std::min(tmp, a[p+w[k]]);
          }
	  b[p] = tmp;
        }
    }
}

void test_dilation_kernel(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));

void
test_dilation_kernel(const image2d<int>& a, image2d<int>& b)
{
  //fill(b, std::numeric_limits<int>::max());
  extension::fill(b, std::numeric_limits<int>::max());

  {
    const kernel::image<0> f {};
    const kernel::image<1> out {};
    using kernel::placeholders::p;
    using kernel::placeholders::n;
    const kernel::Aggregate< mln::accu::features::min<> > Min = {};

    auto expr = (out(p) = Min(f(n)));

    kernel::execute(expr, c8, a, b);
  }
}

void test_extgrad_kernel(const image2d<int>& a, image2d<int>& b) __attribute__ ((noinline));

void
test_extgrad_kernel(const image2d<int>& a, image2d<int>& b)
{
  //fill(b, std::numeric_limits<int>::max());
  extension::fill(b, std::numeric_limits<int>::max());

  {
    const kernel::image<0> f {};
    const kernel::image<1> out {};
    using kernel::placeholders::p;
    using kernel::placeholders::n;

    const kernel::Aggregate< mln::accu::features::min<> > Min = {};
    const kernel::Aggregate< mln::accu::features::max<> > Max = {};

    auto expr = (out(p) = Max(f(n)) - Min(f(n)));

    kernel::execute(expr, c8, a, b);
  }
}



/*
void test_dilation_extfor(const image2d<int>& a, image2d<int>& b)
{
  fill(b, std::numeric_limits<int>::max());

  mln_pixter (px, a);
  mln_pixter (pout, b);
  mln_nixter (nx, c4(px));

  forall(px, pout)
    {
      //std::cout << px.point() << ":" << pout.point() << " {";
      forall(nx)
      {
        //std::cout << nx.val() << ", ";
        pout.val() = std::min(pout.val(), nx.val());
      }
      //std::cout << "}" << std::endl;
    }
}
*/

int main()
{
  int NTEST = 50;

  image2d<int> ima1(1000, 10000);
  image2d<int> ima2(1000, 10000);

  iota(ima1, 0);
  iota(ima2, 1);

  auto zp = imzip(ima1, ima2);

  //io::imprint(imzip(ima1, ima2));

  boost::timer t;
  double r, thistime;

  access_point(ima1, point2d{0,0});
  access_index(ima1, 1000);

  std::cout << "Native pointer..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_native(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;


  std::cout << "Normal Iterators..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_iterator(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;


  std::cout << "Zip Value Iterators..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_zip(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;

  std::cout << "Zip Pixel Iterators..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_zip_pix(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;


  std::cout << "Zip For loop..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_for(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;

  std::cout << "Zip For loop pixel..." << std::endl;
  t.restart();
  r = 0;
  for (int i = 0; i < NTEST; ++i)
    r += test_for_pixel(ima1, ima2);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << " R:" << r <<std::endl;


  std::cout << "Zip For Dilation..." << std::endl;
  image2d<int> ref(1000, 10000);
  t.restart();
  for (int i = 0; i < NTEST/2; ++i)
    test_dilation_pixel(ima1, ref);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime <<std::endl;

  {
    std::cout << "Dilation without zip..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_dilation_pixel_2(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << equal(out, ref) <<std::endl;
  }

  {
    std::cout << "Point Dilation 2..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_dilation_point(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << equal(out, ref) <<std::endl;
  }

  {
    std::cout << "Dilation with index..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_dilation_index(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << equal(out, ref) <<std::endl;
  }

  {
    std::cout << "Native Dilation..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_dilation_native(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << equal(out, ref) <<std::endl;
  }

   {
    std::cout << "Kernel Dilation..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_dilation_kernel(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << equal(out, ref) <<std::endl;
  }


   {
     std::cout << "Kernel Ext-gradient..." << std::endl;
    image2d<int> out(1000, 10000);
    t.restart();
    for (int i = 0; i < NTEST/2; ++i)
      test_extgrad_kernel(ima1, out);
    thistime = t.elapsed();
    std::cout << "Elapsed: " << thistime << " Eq:" << 1 <<std::endl;
  }
}

