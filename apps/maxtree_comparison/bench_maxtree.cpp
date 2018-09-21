#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/range/algorithm/generate.hpp>
#include <mln/morpho/maxtree_ufind.hpp>
#include <mln/morpho/maxtree_ufindrank_parallel.hpp>
#include <mln/morpho/maxtree_hqueue_parallel.hpp>
#include <mln/morpho/maxtree_pqueue_parallel.hpp>
#include <mln/morpho/maxtree_ufind_parallel.hpp>
#include <mln/core/grays.hpp>
#include <mln/io/imprint.hpp>
#include <boost/program_options.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/tick_count.h>

#include <mln/core/value/int.hpp>
#include <random>

// namespace mln
// {


//   void unify_parent(const mln::image2d<uint8>& f,
//                     const std::vector<std::size_t>& S,
//                     mln::image2d<std::size_t>& parent)

//   {
//     mln::image2d< std::set<std::size_t> > ima;
//     resize(ima, parent);

//     mln_foreach(auto& px, parent.pixels())
//       {
//         std::size_t p = px.index();
//         std::size_t q = parent[p];
//         if (f[p] == f[parent[p]]) // p and q belong to the same node.
//           if (p < parent[p])
//             {
//               parent[p] = parent[q];
//               parent[q] = p;
//             }
//       }

//     // Recanonize
//     for (unsigned i = 0; i < S.size(); ++i)
//       {
//         std::size_t p = S[i];
//         std::size_t q = parent[p];
//         if (f[q] == f[parent[q]])
//           parent[p] = parent[q];
//       }
//   }

// }

// bool iscanonized(const mln::image2d<mln::uint8>& ima,
// 		 const mln::image2d<std::size_t>& parent)
// {
//   mln_pixter(px, parent);
//   mln_forall(px)
//   {
//     std::size_t q = px->val();
//     if (not(q == parent[q] or ima[q] != ima[parent[q]]))
//       {
// 	std::cout << "canaonization error @ " << px->index() << std::endl;
// 	return false;
//       }
//   }
//   return true;
// }

// mln::image2d<std::size_t>
// pt2idx(const mln::image2d<mln::point2d>& parent)
// {
//   mln::image2d<std::size_t> out;
//   mln::resize(out, parent);


//   mln_viter(vin, vout, parent, out);
//   mln_forall(vin, vout)
//     *vout = parent.index_of_point(*vin);
//   return out;
// }



// template <typename StrictWeakOrdering>
// void runtest(const mln::image2d<mln::uint8>& ima, StrictWeakOrdering cmp)
// {
//   using namespace mln;

//   image2d<std::size_t> parent1, parent2, parent3, parent4;
//   image2d<point2d> parent4_;
//   std::vector<std::size_t> S1, S2, S3;
//   std::tie(parent1, S1) = morpho::maxtree(ima, c4, cmp);
//   std::tie(parent2, S2) = morpho::maxtree_hqueue(ima, c4, cmp );
//   std::tie(parent3, S3) = morpho::maxtree_ufindbyrank(ima, c4, cmp );
//   parent4_ = morpho::maxtree_ufind_parallel(ima, c4, cmp );
//   parent4 = pt2idx(parent4_);

//   // io::imprint(parent1);
//   // io::imprint(parent4);
//   // io::imprint(parent4_);
//   // io::imprint(parent3);
//   BOOST_CHECK(iscanonized(ima, parent1));
//   BOOST_CHECK(iscanonized(ima, parent2));
//   BOOST_CHECK(iscanonized(ima, parent3));
//   BOOST_CHECK(iscanonized(ima, parent4));

//   unify_parent(ima, S1, parent1);
//   unify_parent(ima, S1, parent2);
//   unify_parent(ima, S3, parent3);
//   unify_parent(ima, S1, parent4);
//   //io::imprint(parent1);
//   //io::imprint(parent3);
//   BOOST_CHECK(all(parent1 == parent2));
//   BOOST_CHECK(all(parent1 == parent3));
//   BOOST_CHECK(all(parent1 == parent4));
// }


int main(int ac, char** av)
{
  using namespace mln;
  using namespace tbb;
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("nthread", po::value<int>()->default_value((int) tbb::task_scheduler_init::automatic), "set number of thread (default: auto)")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }


  typedef UInt<12> V;
  typedef typename image2d<V>::size_type size_type;
  image2d<V> ima(4000, 3500);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> sampler(0, value_traits<V>::max());
  range::generate(ima.values(), [&sampler, &gen] () { return sampler(gen); }) ;

  std::cout << "Number of thread (default): " << tbb::task_scheduler_init::default_num_threads() << std::endl;
  tbb::task_scheduler_init ts(vm["nthread"].as<int>());

  static const int NTEST = 3;

  image2d<size_type> parent;
  image2d<point2d> parent_;
  std::vector<size_type> S;
  std::less<V>     cmp;

  tick_count t0, t1;
  // {
  //   t0 = tick_count::now();
  //   for (int i = 0; i < NTEST; ++i)
  //     std::tie(parent, S) = morpho::maxtree(ima, c4, cmp);
  //   t1 = tick_count::now();
  //   std::cout << "Union Find Maxtree Index: " << (t1-t0).seconds() / NTEST << std::endl;
  // }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::serial::maxtree_ufind(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Serial Union Find Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::parallel::maxtree_ufind(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Parallel Union Find Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::parallel::maxtree_ufind_line(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Parallel Union Find Line Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }


  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::serial::maxtree_hqueue(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Serial Hqueue Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::parallel::maxtree_hqueue(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Parallel Hqueue Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }


  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::serial::maxtree_ufindrank(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Serial Union Find By Rank Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::parallel::maxtree_ufindrank(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Parallel Union Find By Rank Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::serial::maxtree_pqueue(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Serial PQueue Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }

  {
    t0 = tick_count::now();
    for (int i = 0; i < NTEST; ++i)
      std::tie(parent, S) = morpho::impl::parallel::maxtree_pqueue(ima, c4, cmp);
    t1 = tick_count::now();
    std::cout << "Parallel Pqueue Maxtree: " << (t1-t0).seconds() / NTEST << std::endl;
  }



}
