#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/range/algorithm/generate.hpp>
#include <mln/morpho/maxtree_ufind.hpp>
#include <mln/morpho/maxtree_pqueue_2.hpp>

#include <mln/morpho/maxtree_pqueue_parallel.hpp>
#include <mln/morpho/maxtree_ufindrank_parallel.hpp>
#include <mln/morpho/maxtree_hqueue_parallel.hpp>
#include <mln/morpho/maxtree_ufind_parallel.hpp>
#include <mln/morpho/maxtree_najman.hpp>
#include <mln/morpho/maxtree1d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/io/imprint.hpp>
#include <tbb/task_scheduler_init.h>

#define BOOST_TEST_MODULE Morpho
#include <boost/test/unit_test.hpp>
#include <random>

namespace mln
{


  template <typename V, typename size_type>
  void unify_parent(const mln::image2d<V>& f,
                    const std::vector<size_type>&,
                    mln::image2d<size_type>& parent)

  {
    auto px = parent.pixels().riter();
    mln_forall(px)
      {
	//std::cout << px->val() << " ! " << px->index() << std::endl;
        std::size_t p = px->index();
        std::size_t q = mln::morpho::internal::zfind_repr(f, parent, p);
	if (p != q)
	  {
	    if (q == parent[q]) // transmit root property
	      parent[p] = p;
	    else
	      parent[p] = parent[q];
	    parent[q] = p;
	  }
      }

    mln_foreach(auto& p, parent.values())
      p = mln::morpho::internal::zfind_repr(f, parent, p);
  }

}

template <typename V, typename size_type>
bool iscanonized(const mln::image2d<V>& ima,
		 const mln::image2d<size_type>& parent)
{
  mln_pixter(px, parent);
  mln_forall(px)
  {
    std::size_t q = px->val();
    if (not(q == parent[q] or ima[q] != ima[parent[q]]))
      {
	std::cout << "canaonization error @ " << px->index() << std::endl;
	return false;
      }
  }
  return true;
}

template <typename V>
bool iscanonized(const mln::image2d<V>& ima,
		 const mln::image2d<mln::morpho::maxtree_node>& tree)
{
  mln_pixter(px, tree);
  mln_forall(px)
  {
    std::size_t q = px->val().m_parent;
    if (not(q == tree[q].m_parent or ima[q] != ima[tree[q].m_parent]))
      {
	std::cout << "canaonization error @ " << px->index() << std::endl;
	return false;
      }
  }
  return true;
}

bool
check_S(const mln::image2d<mln::morpho::maxtree_node>& tree, unsigned root)
{
  using namespace mln;
  image2d<bool> dejavu;
  resize(dejavu, tree).init(false);

  dejavu[root] = true;
  for (unsigned x = root; x != (unsigned) -1; x = tree[x].m_next)
    {
      assert(dejavu[tree[x].m_parent]);
      if (!dejavu[tree[x].m_parent])
	return false;
      dejavu[x] = true;
    }
  return true;
}



mln::image2d<std::size_t>
pt2idx(const mln::image2d<mln::point2d>& parent)
{
  mln::image2d<std::size_t> out;
  mln::resize(out, parent);


  mln_viter(vin, vout, parent, out);
  mln_forall(vin, vout)
    *vout = parent.index_of_point(*vin);
  return out;
}



template <typename V, typename StrictWeakOrdering>
void runtest(const mln::image2d<V>& ima, StrictWeakOrdering cmp)
{
  using namespace mln;

  typedef typename image2d<V>::size_type size_type;
  image2d<size_type> parent1, parent;
  std::vector<size_type> S1, S;
  std::tie(parent1, S1) = morpho::impl::serial::maxtree_ufind(ima, c4, cmp);
  unify_parent(ima, S1, parent1);

  {
    std::tie(parent, S) = morpho::impl::serial::maxtree_hqueue(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }


  {
    std::tie(parent, S) = morpho::maxtree_najman(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    std::tie(parent, S) = morpho::impl::parallel::maxtree_hqueue(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    std::tie(parent, S) = morpho::impl::serial::maxtree_ufind(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    std::tie(parent, S) = morpho::impl::parallel::maxtree_ufind(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    std::tie(parent, S) = morpho::impl::parallel::maxtree_ufind_line(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }


  {
    std::tie(parent, S) = morpho::impl::serial::maxtree_ufindrank(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    std::tie(parent, S) = morpho::impl::parallel::maxtree_ufindrank(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }



  {
    std::tie(parent, S) = morpho::impl::serial::maxtree_pqueue(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

  {
    image2d<morpho::maxtree_node> tree;
    unsigned root;
    std::tie(tree, root) = morpho::maxtree_pqueue_2(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, tree));
    BOOST_CHECK(check_S(tree, root));
    auto parent = transform(tree, std::mem_fn(&morpho::maxtree_node::m_parent));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }


  {
    std::tie(parent, S) = morpho::impl::parallel::maxtree_pqueue(ima, c4, cmp );
    BOOST_CHECK(iscanonized(ima, parent));
    BOOST_CHECK(morpho::check_S(parent, S.data(), S.data() + S.size()));
    unify_parent(ima, S1, parent);
    BOOST_CHECK(all(parent == parent1));
  }

}


BOOST_AUTO_TEST_CASE(Maxtree)
{
  using namespace mln;
  typedef UInt<8> V;
  typedef typename image2d<V>::size_type size_type;
  image2d<V> ima(300, 100);


  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> sampler(0, value_traits<V>::max());
  range::generate(ima.values(), [&sampler, &gen] () { return sampler(gen); }) ;

  tbb::task_scheduler_init ts;

  // {
  //   image2d<size_type> f;
  //   resize(f, ima);
  //   mln_foreach(auto px, f.pixels())
  //     px.val() = px.index();
  //   io::imprint(f);
  // }

  // io::imprint(ima);

  // {
  //   image2d<point2d> parent_;
  //   image2d<size_type> parent1;
  //   std::vector<size_type> S;
  //   std::less<uint8> cmp;
  //   resize(parent_, ima);

  //   morpho::maxtree1d(ima, parent_, 0, cmp);
  //   std::tie(parent1, S) = morpho::maxtree(ima, c4, cmp);
  //   auto parent = pt2idx(parent_);

  //   BOOST_CHECK(iscanonized(ima, parent));
  //   BOOST_CHECK(iscanonized(ima, parent1));

  //   unify_parent(ima, S, parent1);
  //   unify_parent(ima, S, parent);
  //   BOOST_CHECK(all(parent == parent1));
  // }

  runtest(ima, std::less<V> ());
  runtest(ima, std::greater<V> ());

}
