#include <mln/core/image/image2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/wrt_offset.hpp>
#include <algorithm>
#include <numeric>
#include <boost/timer.hpp>


using namespace mln;


double test_viter(const image2d<int>& ima)
{
  double v = 0;
  mln_foreach(auto& x, ima.values()) {
    v += x;
  }
  return v;
}

double test_pixter(const image2d<int>& ima)
{
  double v = 0;
  mln_foreach(auto& x, ima.pixels())
    v += x.val();

  return v;
}

double test_piter(const image2d<int>& ima)
{
  double v = 0;
  mln_foreach(auto p, ima.domain()) {
    v += ima(p);
  }
  return v;
}



double test_native(const int* ima2, int nrows, int ncols)
{
  double v = 0;
  int end = nrows * ncols;
  for (int i = 0; i < end; ++i)
    v += ima2[i];

  return v;
}


double test_nbh_pixter(const image2d<int>& ima)
{
  double u = 0;
  image2d<int>::const_pixel_range::iterator px = ima.pixels().iter();
  auto nx = c8(*px).iter();


  mln_forall(px)
    mln_forall(nx)
      u += nx->val();

  return u;
}

double test_nbh_piter(const image2d<int>& ima)
{
  double u = 0;
  auto p = ima.domain().iter();
  auto n = c8(*p).iter();

  mln_forall(p)
    mln_forall(n)
	u += ima.at(*n);

  return u;
}

double test_nbh_index(const image2d<int>& ima)
{
  double u = 0;
  std::size_t idx = ima.index_of_point(ima.domain().pmin);
  auto w = wrt_delta_index(ima, c8_t::dpoints);

  //std::cout << ima.index_strides()[0] << "," << ima.index_strides()[1] << std::endl;

  unsigned nrows = ima.nrows();
  unsigned ncols = ima.ncols();
  for (unsigned i = 0; i < nrows; ++i)
    {
      for (unsigned j = 0; j < ncols; ++j)
        {
          std::size_t p = idx + j;// * ima.index_strides()[1];
          mln_foreach(auto k, w) {
            u += ima[p + k];
          }
        }
      idx += ima.index_strides()[0];
    }

  return u;
}


double test_native_nbh(const image2d<int>& ima)
{
  double r2 = 0;

  const int sz = 8;
  auto dpoints = c8_t::dpoints;


  const char* ptr2 = (char*) &ima(ima.domain().pmin);
  const size_t* strides = ima.strides();
  const point2d pmax = ima.domain().shape();
  ptrdiff_t offsets[8];
  wrt_offset(ima, dpoints, offsets);

  for (int x = 0; x < pmax[0]; ++x, ptr2 += strides[0])
    {
      const char* ptr = ptr2;
      for (int y = 0; y < pmax[1]; ++y, ptr += strides[1])
        for (int k = 0; k < sz; ++k)
          r2 += *(const int*)(ptr + offsets[k]);
    }

  return r2;
}


template <typename T>
void iota(image2d<T>& ima, T v)
{
  mln_foreach(auto& x, ima.values())
    x = v++;
}


void display()
{
  const int nrows = 5, ncols = 5;
  image2d<int> ima(nrows, ncols);
  iota(ima, 0);

  {
    std::cout << "Display forward site iterator." << std::endl;
    mln_foreach(auto p, ima.domain())
      std::cout << p << ",";
    std::cout << std::endl;
  }
  {
    std::cout << "Display forward value iterator." << std::endl;
    mln_viter(x, ima);
    mln_forall(x)
      std::cout << *x << ",";
    std::cout << std::endl;
  }
  {
    std::cout << "Display forward pixel iterator." << std::endl;
    mln_pixter(x, ima);
    mln_forall(x)
      std::cout << "(" << x->point() << ":" << x->val() << "," << x->index() << "),";
    std::cout << std::endl;
  }
}

void display_nbh()
{
  const int nrows = 5, ncols = 5;
  image2d<int> ima(nrows, ncols);
  iota(ima, 0);

  {
    std::cout << "Display forward site iterator." << std::endl;
    auto p = ima.domain().iter();
    mln_iter(n, c8(*p));
    mln_forall(p)
    {
      std::cout << *p << ": ";
      mln_forall(n)
        std::cout << *n << ",";
      std::cout << std::endl;
    }
  }

  {
    std::cout << "Display forward pixel iterator." << std::endl;
    mln_pixter(px, ima);
    mln_iter(nx, c8(*px));

    mln_forall(px)
    {
      std::cout << "{" << px->point() << "," << px->val() << "," << px->index() << "}: ";
      mln_forall(nx)
        std::cout << "{" << nx->point() << "," << nx->val() << "}: ";
      std::cout << std::endl;
    }
  }
}


int main()
{

  const int nrows = 1000, ncols = 10000;
  image2d<int> ima(nrows, ncols);
  iota(ima, 0);

  display();
  display_nbh();

  //std::iota(std::begin(ima.values()), std::end(ima.values()), 1);

  static const int ntest = 10;

  boost::timer t;
  double r, r1, r2, r3, thistime;

  std::cout << "Point Iterators..." << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r  = test_piter(ima);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << std::endl;
  std::cout << "Result: " << r << std::endl;

  std::cout << "Value Iterators..." << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r1  = test_viter(ima);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << std::endl;
  std::cout << "Result: " << r1 << std::endl;


  std::cout << "Pixel Iterators..." << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r3  = test_pixter(ima);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << std::endl;
  std::cout << "Result: " << r3 << std::endl;

  // defaut with pointers
  int* ima2 = new int[nrows * ncols];
  std::iota(ima2, ima2 + nrows * ncols, 1);


  std::cout << "Native..." << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r2 = test_native(ima2, nrows, ncols);
  thistime = t.elapsed();
  std::cout << "Elapsed: " << thistime << std::endl;

  std::cout << r1 << " " << r2 << " " << r3 << std::endl;


  // bench neighborhood


  std::cout << "Neighborhood... piter/niter" << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r1 = test_nbh_piter(ima);
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
  std::cout << r1 << std::endl;


  std::cout << "Neighborhood... pixter iterators" << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r1 = test_nbh_pixter(ima);
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
  std::cout << r1 << std::endl;


  std::cout << "Neighborhood... indexes" << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r1 = test_nbh_index(ima);
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
  std::cout << r1 << std::endl;


  std::cout << "Neighborhood native..." << std::endl;
  t.restart();
  for (int i = 0; i < ntest; ++i)
    r1 = test_native_nbh(ima);
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
  std::cout << r1 << std::endl;

  delete [] ima2;
}


