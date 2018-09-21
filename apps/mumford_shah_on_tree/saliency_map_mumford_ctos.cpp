#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/morpho/component_tree/io.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/extinction.hpp>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/dontcare.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/property_map/function_property_map.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/croutines.hpp>
#include <apps/g2/compute_g2.hpp>
#include <apps/g2/satmaxtree.hpp>
#include "mumford_shah_on_tree.hpp"

// Compute the depth attribute of each graph node
boost::vector_property_map<unsigned>
compute_graph_depth(const MyGraph& g)
{
  mln_entering("Compute graph depth");

  boost::vector_property_map<unsigned> depth(boost::num_vertices(g));

  auto one = [](mln::dontcare_t) -> int{ return 1; };
  auto w = boost::make_function_property_map<MyGraph::edge_descriptor, int, decltype(one)>(one);

  MyGraph::vertex_descriptor root = boost::vertex(0, g);
  depth[root] = 0;

  MyGraph gT;
  boost::transpose_graph(g, gT);
  boost::dag_shortest_paths(gT, root, boost::weight_map(w)
                            .distance_map(depth)
                            .distance_compare(std::greater<int> ())
                            .distance_inf(-1)
                            .distance_zero(0)
                            );
  mln_exiting();
  return depth;
}

// Compute the per-pixel attribute and reconstruct
template <class ValueMap>
void
write_vmap_to_image(const tree_t* t, const tlink_t* tlink,
                    const ValueMap& vmap, mln::image2d<mln::uint16>& out)
{
  mln_foreach(auto px, out.pixels())
  {
    unsigned w = 0;
    for (int i = 0; i < NTREE; ++i)
      {
        tree_t::node_type tnode = t[i].get_node_at(px.index());
        MyGraph::vertex_descriptor gnode = tlink[i][tnode];
        w = std::max(w, vmap[gnode]);
      }
    px.val() = w;
  }
}

/// \brief Remove non-2F from the tree
template <class P>
mln::morpho::component_tree<P, mln::image2d<P> >
tree_keep_2F(const mln::morpho::component_tree<P, mln::image2d<P> >& tree)
{
  using namespace mln;
  morpho::component_tree<P, image2d<P> > out;

  auto newdata = out._get_data();
  auto olddata = tree._get_data();

  // 1. Copy the point2node map
  box2d olddom = olddata->m_pmap.domain();
  box2d dom;
  dom.pmin = olddom.pmin / 2;
  dom.pmax = (olddom.pmax + 1) / 2;
  newdata->m_pmap.resize(dom);
  copy(olddata->m_pmap | sbox2d{olddom.pmin, olddom.pmax, {2,2}},
       newdata->m_pmap);

  // 2. Copy the node
  newdata->m_nodes = olddata->m_nodes;

  // 3. Copy the point set and update node first point index/
  newdata->m_S.resize(dom.size());
  unsigned j = 0;
  for (unsigned i = 0; i < olddata->m_S.size(); ++i)
    {
      P p = olddata->m_S[i];
      point2d pt = olddata->m_pmap.point_at_index(p);
      if (K1::is_face_2(pt))
        {
          newdata->m_S[j] = newdata->m_pmap.index_of_point(pt/2);
          auto node = tree.get_node_at(p);
          if (node.get_first_point_id() == i)
            newdata->m_nodes[node.id()].m_point_index = j;
          ++j;
        }
    }
  // 4. Do not forget the sentinel
  newdata->m_nodes[out.npos()].m_point_index = j;

  return out.get_subtree(tree.get_root_id());
}



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] α₀ α₁ output[float]\n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n";
  std::exit(1);
}


int main(int argc, char** argv)
{
  if (argc < 4)
    usage(argv);



  const char* input_path = argv[1];
  int a0 = std::atoi(argv[2]);
  int a1 = std::atoi(argv[3]);
  const char* output_path = argv[4];


  tbb::task_scheduler_init init;


  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<V> ima;
  io::imread(input_path, ima);

  image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
  image2d<V> F = interpolate_k1(f);


  // COmpute the color tree of shapes T
  tree_t T;
  {
    // 1. Compute the marginal ToS and filter them if necessary.
    tree_t t[NTREE];
    tbb::parallel_for(0, (int)NTREE, [&t,&f,a0](int i){
        t[i] = morpho::cToS(imtransform(f, [i](value_t x) { return x[i]; }), c4);
        if (a0 > 0) {
          grain_filter_inplace(t[i], a0);
          t[i].shrink_to_fit();
        }
      });

    // 2. Compute the Gos.
    MyGraph g2;
    std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
    std::tie(g2, tlink) = compute_g2<NTREE>(t);

    // 3. Compute the depth image
    boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
    image2d<uint16> imdepth = imchvalue<uint16>(F);

    write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

    // debug
    io::imsave(imdepth, "depth.tiff");

    // 4. Compute the saturated maxtree
    std::tie(T, std::ignore) = satmaxtree(imdepth);
    T = tree_keep_2F(T);
  }

  //morpho::save(T, "out.tree");

  // Mumford-shah + simplification
  {
    // Grain filter
    if (a1 > 0)
      {
        grain_filter_inplace(T, a1);
        T.shrink_to_fit();
      }

    //std::cout << "computation of saliency " << std::endl;
    auto aminnu = compute_nu(T, F);

    // mln_foreach(tree_t::node_type x, T.nodes())
    //   aenergy[x] = aminnu[x];

    auto ienergy = morpho::make_image(T, aminnu);

    auto extincted = morpho::extinction(ienergy, morpho::tree_neighb_t());
    property_map<tree_t, float> vmap = std::move(extincted.get_vmap());
    image2d<float> sal = set_value_on_contour(T, vmap);
    //io::imsave(sal, "test.tiff");
    io::imsave(sal, output_path);
    // std::cout << "Simplification (number of llines): " << std::endl;
    // mumford_shah_on_tree(T, F, lambda);
    // F = Kadjust_to(F, ima.domain());
    // io::imsave(F, output_path);
  }


}
