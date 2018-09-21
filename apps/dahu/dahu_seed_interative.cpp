
#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>

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



#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"


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



namespace mln
{

    bool is_top (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 )
            return true;
        else
            return false;

    }

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] != 0  and p[0]  != 1  and p[0] != height -1 and p[0] != height -2 and p[1] != 0  and p[1]  != 1 and p[1]  != width - 1 and p[1] != width -2)
        {
            if (p[0] == 2 or p[0] == height - 3 or p[1] == width -3 or p[1] == 2  )
                return true;
            else
                return false;
        }

        //if (p[0] == 0 or p[0] == height - 1 or p[1] == width -1 or p[1] == 0  )
        //	return true;
        //else
        //	return false;

    }
}



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] input1[rgb] α₀ α₁ λ output[rgb]\n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}


int main(int argc, char** argv)
{
  if (argc < 6)
    usage(argv);



  const char* input_path = argv[1];
  const char* seed_path = argv[2];

  int a0 = std::atoi(argv[3]);
  int a1 = std::atoi(argv[4]);
  int lambda = std::atoi(argv[5]);
  const char* output_path = argv[6];


  tbb::task_scheduler_init init;





  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<V> ima;
  io::imread(input_path, ima);

  image2d<V> image_seed;
  io::imread(seed_path,image_seed);
  image2d<uint8_t> image_seed_gray(image_seed.domain());


  mln_foreach(point2d p1, image_seed.domain())
  {
      image_seed_gray(p1) = 0.2126 * image_seed(p1)[0]+ 0.7152* image_seed(p1)[1] + 0.0722 * image_seed(p1)[2];
  }
  image2d<uint8_t> seed_image = interpolate_k1(image_seed_gray);


  image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
  image2d<V> F = interpolate_k1(f);
  unsigned height = F.nrows();
  unsigned width = F.ncols();

  //std::string outputDirectory = "/lrde/home/movn/Documents/code/code_edwin2/program/apps/mumford_shah_on_tree/abc";

  image2d<V> ima_compo(F.domain() );



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


    auto& U0  = t[0]._get_data()->m_Uv;
    auto& U1  = t[1]._get_data()->m_Uv;
    auto& U2  = t[2]._get_data()->m_Uv;


    mln_foreach(point2d p1, U0.domain())
    {
        ima_compo(p1)[0] = U0(p1);
        ima_compo(p1)[1] = U1(p1);
        ima_compo(p1)[2] = U2(p1);
    }

    // 2. Compute the Gos.
    MyGraph g2;
    std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
    std::tie(g2, tlink) = compute_g2<NTREE>(t);

    // 3. Compute the depth image
    boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
    image2d<uint16> imdepth = imchvalue<uint16>(F);

    write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

    // debug
    // io::imsave(imdepth, "depth.tiff");

    // 4. Compute the saturated maxtree
    std::tie(T, std::ignore) = satmaxtree(imdepth);
    T = tree_keep_2F(T);
    T._reorder_pset();





    std::ofstream file("abc.dot",
             std::ios_base::out | std::ios_base::trunc);
    write_graphviz(file, T);

  }




    unsigned nodenumbers = T.nodes().size();

    auto& nodes1 = T._get_data()->m_nodes;
    auto& S     = T._get_data()->m_S;
    auto& pmap  = T._get_data()->m_pmap;

    box2d D = pmap.domain();

    image2d <uint8_t>  depth(D);
    extension::fill(depth, 0);







    // dahu_distance




    // initiation


    std::vector<V>  min_node(nodenumbers);
    std::vector<V>  max_node(nodenumbers);

    std::vector<V>  min_node_final(nodenumbers);
    std::vector<V>  max_node_final(nodenumbers);

    std::vector<V>  dist(nodenumbers);

    image2d<V> distancemap(D );
    image2d<V> min_pixel(D );
    image2d<V> max_pixel(D );
    image2d<uint8_t> distancemap_final_gray(D );
    image2d<uint8_t> distancemap_final_gray1(D );


    std::vector< std::vector<unsigned> > children(nodenumbers);


    //int dejavu_node[nodenumbers] = {0};
    //dejavu_node[pmap(seedpoints)]  = 1;


    //int seed_node[nodenumbers] = {0};
    std::vector<int>  seed_node(nodenumbers);
    std::fill (seed_node.begin(),seed_node.end(),0);

    seed_node[0]  = 1;

    std::vector<int>  dejavu_node(nodenumbers);
    std::fill (dejavu_node.begin(),dejavu_node.end(),0);

    //int dejavu_node[nodenumbers] = {0};
    dejavu_node[0]  = 1;

    image2d<uint8_t> test(D );

    mln_foreach(auto p, D)
    {
        if (seed_image(p) >250)
        {
            dejavu_node[pmap(p)] =1;
            seed_node[pmap(p)] =1;
            test(p) = 255;
        }
    }


    //std::cout << "seed node   "  << pmap(seedpoints)  << std::endl;




    std::cout << "initial  step"   << std::endl;

    mln_foreach(auto x, T.nodes())
    {


        std::vector<V> listpixel1;

        mln_foreach (auto p, x.proper_pset())
        {
            listpixel1.push_back(ima_compo(ima_compo.point_at_index(p)));

        }


        std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end(),lexicographicalorder_less<value_t>());

        V medianvalue1;

        medianvalue1 = listpixel1[listpixel1.size()/2];
        //std::cout << "median value"  << medianvalue << "  " << medianvalue1  << std::endl;


        for (int i = 0; i< 3; ++i)
        {
            min_node[x.id()][i] = medianvalue1[i];
            max_node[x.id()][i] = medianvalue1[i];

            min_node_final[x.id()][i] = medianvalue1[i];
            max_node_final[x.id()][i] = medianvalue1[i];
            if (seed_node[x.id() == 1])
                dist[x.id()][i] = 0;
            else
                dist[x.id()][i] = 255;
        }

        unsigned q = x.get_parent_id();
        children[q].push_back(x.id());

    }





    // test value
    /*
    std::cout << "test value  "  << std::endl;
    mln_foreach(auto x, T.nodes())
    {
        std::cout << "node "  << x.id()  << std::endl;
        for (int i = 0; i < 3; ++i)
        {
            std::cout << "value  " << int(min_node[x.id()][i])  << std::endl;
        }
    }

    */


    std::cout << "chay up len " << std::endl;

    mln_reverse_foreach(auto x, T.nodes())
    {
        if (seed_node[x.id()] == 0)
        {

            const unsigned nchildren = children[x.id()].size();

            std::vector <unsigned > dist_temp(nchildren);
            std::fill (dist_temp.begin(),dist_temp.end(),0);

            std::vector <V> min_temp(nchildren);
            std::vector <V> max_temp(nchildren);

            std::fill (min_temp.begin(),min_temp.end(),min_node[x.id()]);
            std::fill  (max_temp.begin(),max_temp.end(),max_node[x.id()]);

            int min_dis = 255 * 3;

            for (int j = 0 ; j < nchildren ; ++j)
            {
                unsigned child = children[x.id()][j];

                if (dejavu_node[child] == 1)
                {
                    for (int i = 0 ; i < 3 ; ++i)
                    {
                        if (min_temp[j][i]  > min_node[child][i])
                            min_temp[j][i] = min_node[child][i];
                        if (max_temp[j][i] < max_node[child][i])
                            max_temp[j][i] = max_node[child][i];

                        dist_temp[j] = dist_temp[j] + max_temp[j][i] - min_temp[j][i];

                    }
                    dejavu_node[x.id()] = 1;
                }
                else
                    dist_temp[j] = 255 *3 ;

                if (dist_temp[j] < min_dis)
                {
                    min_dis = dist_temp[j];
                    min_node[x.id()] = min_temp[j];
                    max_node[x.id()] = max_temp[j];
                }
            }
        }
    }



    std::cout << "chay down xuong"  << std::endl;

    mln_foreach(auto x, T.nodes())
    {

        if (seed_node[x.id()] == 0)
        {
            V min_temp = min_node_final[x.id()];
            V max_temp = max_node_final[x.id()];
            if (dejavu_node[x.id()] == 0)
            {
                for (int i = 0; i < 3 ; ++i)
                {
                    if (min_node[x.id()][i] > min_node[x.get_parent_id()][i] )
                        min_node[x.id()][i] = min_node[x.get_parent_id()][i];
                    if (max_node[x.id()][i] < max_node[x.get_parent_id()][i])
                        max_node[x.id()][i] = max_node[x.get_parent_id()][i];
                }
                dejavu_node[x.id()] = 1;
            }
            else
            {
                int dist_temp = 0;
                int dist_prev = 0;
                for (int i = 0 ; i < 3 ; ++i)
                {
                    if (min_temp[i] > min_node[x.get_parent_id()][i] )
                        min_temp[i] = min_node[x.get_parent_id()][i];
                    if (max_temp[i] < max_node[x.get_parent_id()][i])
                        max_temp[i] = max_node[x.get_parent_id()][i];
                    dist_temp = dist_temp + max_temp[i] - min_temp[i];
                    dist_prev = dist_prev + max_node[x.id()][i] - min_node[x.id()][i];
                }
                if (dist_temp < dist_prev)
                {
                    min_node[x.id()] = min_temp;
                    max_node[x.id()] = max_temp;
                }
            }
        }
    }

    mln_foreach(auto x, T.nodes())
    {
        if (seed_node[x.id()] == 1)
            std::cout << max_node[x.id()]  << "   "  << min_node[x.id()]  << std::endl;
    }

    //  back propagation

    mln_foreach(auto x, T.nodes())
    {
        mln_foreach (auto p, x.proper_pset())
        {

            min_pixel(F.point_at_index(p)) = min_node[x.id()];
            max_pixel(F.point_at_index(p)) = max_node[x.id()];
            distancemap(F.point_at_index(p)) = 	max_pixel(F.point_at_index(p)) -  min_pixel(F.point_at_index(p)) ;
            distancemap_final_gray(F.point_at_index(p)) =  distancemap(F.point_at_index(p))[0]/3+ 	distancemap(F.point_at_index(p))[1]/3 + distancemap(F.point_at_index(p))[2]/3		;

        }
    }


    io::imsave(distancemap_final_gray, output_path);

    io::imsave(test, "vkl.pgm");









}
