#include <mln/core/domain/box.hpp>
#include <mln/core/ndimage.hpp>
#include <boost/timer.hpp>
#include <boost/range/algorithm_ext/iota.hpp>
#include <mln/core/grays.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

void bench_home(mln::box3d b)
{
  boost::timer t;

  std::cout << ">>> home made" << std::endl;
  mln::point3d r;
  t.restart();

  for (int i = 0; i < 10; ++i)
    r = std::accumulate(b.rbegin(), b.rend(), mln::point3d{0,0,0});
  std::cout << "Res: " << r << std::endl;
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
}

template <typename T>
void bench_home_ima(const mln::image3d<T>& ima)
{
  boost::timer t;

  std::cout << ">>> home made" << std::endl;
  double r;
  t.restart();
  auto rbegin = boost::begin(ima.rvalues());
  auto rend = boost::end(ima.rvalues());
  for (int i = 0; i < 10; ++i)
    r = std::accumulate(rbegin, rend, 0.0);
  std::cout << "Res: " << r << std::endl;
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
}


void bench_spe(mln::box3d b)
  {
    boost::timer t;

    std::cout << ">>> std::reverse" << std::endl;
    mln::point3d r;
    t.restart();
    boost::reverse_iterator< mln::box3d::iterator > rbegin(b.end()), rend(b.begin());
    for (int i = 0; i < 10; ++i)
      r = std::accumulate(rbegin, rend, mln::point3d{0,0,0});
    std::cout << "Res: " << r << std::endl;
    std::cout << "Elapsed: " << t.elapsed() << std::endl;
  }

template <typename T>
void bench_boost_ima(const mln::image3d<T>& ima)
{
  boost::timer t;

  std::cout << ">>> boost::reverse" << std::endl;
  double r;
  t.restart();
  boost::reverse_iterator< typename mln::image3d<T>::const_value_iterator >
    rbegin(boost::end(ima.values())),
    rend(boost::begin(ima.values()));
  for (int i = 0; i < 10; ++i)
    r = std::accumulate(rbegin, rend, 0.0);
  std::cout << "Res: " << r << std::endl;
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
}

template <typename T>
void bench_fwd_ima(const mln::image3d<T>& ima)
{
  boost::timer t;

  std::cout << ">>> FWD" << std::endl;
  double r;
  t.restart();
  for (int i = 0; i < 10; ++i)
    r = std::accumulate(boost::begin(ima.values()), boost::end(ima.values()), 0.0);
  std::cout << "Res: " << r << std::endl;
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
}


void bench_std(mln::box3d b)
  {
    boost::timer t;

    std::cout << ">>> std::reverse" << std::endl;
    mln::point3d r;
    t.restart();
    std::reverse_iterator< mln::box3d::iterator > rbegin(b.end()), rend(b.begin());
    for (int i = 0; i < 10; ++i)
      r = std::accumulate(rbegin, rend, mln::point3d{0,0,0});
    std::cout << "Res: " << r << std::endl;
    std::cout << "Elapsed: " << t.elapsed() << std::endl;
  }

template <typename T>
void bench_std_ima(const mln::image3d<T>& ima)
{
  boost::timer t;

  std::cout << ">>> std::reverse" << std::endl;
  double r;
  t.restart();
  std::reverse_iterator< typename mln::image3d<T>::const_value_iterator >
    rbegin(boost::end(ima.values())),
    rend(boost::begin(ima.values()));
  for (int i = 0; i < 10; ++i)
    r = std::accumulate(rbegin, rend, 0.0);
  std::cout << "Res: " << r << std::endl;
  std::cout << "Elapsed: " << t.elapsed() << std::endl;
}



int main()
{

  mln::box3d b = {{3,18,25}, {3000, 450, 79}};
  std::cout << "== std::reverse vs home made rev_iterator ==" << std::endl;

  bench_home(b);
  bench_spe(b);
  bench_std(b);


  mln::image3d<mln::uint8> ima(b);
  boost::iota(ima.values(), 0);

  bench_fwd_ima(ima);
  bench_home_ima(ima);
  bench_boost_ima(ima);
  bench_std_ima(ima);

}
