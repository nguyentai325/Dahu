#ifndef MUMFORD_SHAH_ON_TREE_HPP
# define MUMFORD_SHAH_ON_TREE_HPP

# include <mln/morpho/component_tree/accumulate.hpp>
# include <mln/morpho/component_tree/reconstruction.hpp>
# include <mln/morpho/component_tree/filtering.hpp>
# include <mln/accu/accumulators/variance.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/accu/accumulators/accu_as_it.hpp>
# include <apps/attributes/gradient_magnitude.hpp>
# include <mln/core/trace.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>

/// \brief Perform the MS simplification
///
/// \param[inout] tree The input tree
/// \param[inout] F The input image (interpolated to match tree's domain)
/// \param lambda Mumford-shash regularization parameter
template <class T, class I>
void
mumford_shah_on_tree(T& tree,
                     I& F,
                     double lambda);

template <class T, class I>
mln::property_map<T, float>
compute_nu(T& tree, I& F);

/***************************************/
/*** Implementation                 ****/
/***************************************/

namespace internal
{
  typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

  tree_t::vertex_id_t
  find_canonical(mln::property_map<tree_t, tree_t::vertex_id_t>& newpar,
                 tree_t::vertex_id_t x)
  {
    if (newpar[x] == tree_t::npos()) // x is alive
      return x;
    else
      return newpar[x] = find_canonical(newpar, newpar[x]);
  }

}


template <class T, class I>
void
mumford_shah_on_tree(T& tree,
                     I& F,
                     double lambda)
{
  using namespace mln;
  typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

  auto A = morpho::vaccumulate_proper(tree, F,
                                      accu::accumulators::accu_as_it<
                                        accu::accumulators::variance< rgb8, rgb<double>, rgb<double> >
                                      > ());
  auto B = compute_attribute_on_contour(tree, F, accu::features::count<>());

  // Sort the node by gradient
  auto mgrad = compute_gradient_magnitude(tree, F); // note: minima = strong gradient
  std::vector<tree_t::vertex_id_t> S;
  S.reserve(tree.size());

  mln_foreach(tree_t::node_type x, tree.nodes_without_root())
    S.push_back(x.id());

  std::sort(S.begin(), S.end(), [&mgrad](tree_t::vertex_id_t x, tree_t::vertex_id_t y) {
      return mgrad[x] > mgrad[y];
    });


  // Energy function
  auto denergy = [lambda, &A, &B](tree_t::node_type x, tree_t::node_type q) {
    auto tmp = A[q];
    double before = (tmp.to_result() * accu::extractor::count(tmp));

    tmp.take(A[x]);
    double after = (tmp.to_result() * accu::extractor::count(tmp));
    //std::cout << "before: " << before << " after: " << after << std::endl;
    double delta = -lambda * B[x] + (after - before);
    //std::cout << "delta: " << delta << std::endl;
    return delta;
  };

  // Minimization glutone
  std::cout << "Before: " << tree.size() << std::endl;

  unsigned k = tree.size();
  property_map<tree_t, bool> alive(tree, true);
  property_map<tree_t, tree_t::vertex_id_t> newpar(tree, tree_t::npos());

  for (tree_t::vertex_id_t i : S)
    {
      tree_t::node_type x = tree.get_node(i);
      tree_t::node_type q = tree.get_node(::internal::find_canonical(newpar, x.get_parent_id()));

      if (denergy(x,q) < 0)
        {
          alive[x] = false;
          A[q].take(A[x]);
          newpar[x] = q.id();
          --k;
        }
    }

  std::cout << "After: " << k << std::endl;


  {
    morpho::filter_direct_inplace(tree, alive);
    tree.shrink_to_fit();

    auto vmap = morpho::vaccumulate_proper(tree, F, accu::features::mean<> ());
    morpho::reconstruction(tree, vmap, F);
  }

}




template <class T, class I>
mln::property_map<T, float>
compute_nu(T& tree, I& F)
{
  using namespace mln;
  typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

  auto A = morpho::vaccumulate_proper(tree, F,
                                      accu::accumulators::accu_as_it<
                                        accu::accumulators::variance< rgb8, rgb<double>, rgb<double> >
                                      > ());
  auto B = compute_attribute_on_contour(tree, F, accu::features::count<>());

  // Sort the node by gradient
  auto mgrad = compute_gradient_magnitude(tree, F); // note: minima = strong gradient
  std::vector<tree_t::vertex_id_t> S;
  S.reserve(tree.size());

  mln_foreach(tree_t::node_type x, tree.nodes_without_root())
    S.push_back(x.id());

  std::sort(S.begin(), S.end(), [&mgrad](tree_t::vertex_id_t x, tree_t::vertex_id_t y) {
      return mgrad[x] > mgrad[y];
    });


  // Energy function
  auto min_nu = [&A, &B](tree_t::node_type x, tree_t::node_type q) {
    auto tmp = A[q];
    double before = (tmp.to_result() * accu::extractor::count(tmp));

    tmp.take(A[x]);
    double after = (tmp.to_result() * accu::extractor::count(tmp));
    //std::cout << "before: " << before << " after: " << after << std::endl;
    //double delta = -lambda * B[x] + (after - before);
    double nu = (after - before)/B[x];
    //std::cout << "delta: " << delta << std::endl;
    return nu;
  };

  // Minimization glutone
  //std::cout << "Before: " << tree.size() << std::endl;

  property_map<tree_t, float> out(tree);

  unsigned k = tree.size();
  property_map<tree_t, bool> alive(tree, true);
  property_map<tree_t, tree_t::vertex_id_t> newpar(tree, tree_t::npos());

  for (tree_t::vertex_id_t i : S)
  {
    tree_t::node_type x = tree.get_node(i);
    tree_t::node_type q = tree.get_node(::internal::find_canonical(newpar, x.get_parent_id()));

    out[x] = min_nu(x,q);
  }

  for (tree_t::vertex_id_t i : S)
  {
    tree_t::node_type x = tree.get_node(i);
    tree_t::node_type q = tree.get_node(::internal::find_canonical(newpar, x.get_parent_id()));

    out[x] = std::max(out[x], (float)min_nu(x,q));

    //if (denergy(x,q) < 0)
    //{
    alive[x] = false;
    A[q].take(A[x]);
    newpar[x] = q.id();
    --k;
    //}
  }

  float max_nu = 0;
  for (tree_t::vertex_id_t i : S)
  {
    tree_t::node_type x = tree.get_node(i);
    max_nu = std::max(max_nu, out[x]);
  }

  for (tree_t::vertex_id_t i : S)
  {
    tree_t::node_type x = tree.get_node(i);
    out[x] = max_nu - out[x];
    //out[x] = std::log10(1 + out[x]);
  }
  //std::cout << "max_nu = " << max_nu << std::endl;
  //std::cout << "After: " << k << std::endl;


  // {
  //   morpho::filter_direct_inplace(tree, alive);
  //   tree.shrink_to_fit();

  //   auto vmap = morpho::vaccumulate_proper(tree, F, accu::features::mean<> ());
  //   morpho::reconstruction(tree, vmap, F);
  // }



  //***********
  //property_map<tree_t, float> out(tree);
  // mln_foreach(auto x, tree.nodes())
  //   out[x] = f(x);

  // mln_exiting();
  return out;
  //***********

}


#endif // ! MUMFORD_SHAH_ON_TREE_HPP
