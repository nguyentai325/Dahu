
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

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 or  p[0] == height -1 or p[1] == 0  or p[1] == width -1 )
            return true;
        else
            return false;
    }

    template <typename V>
    image2d<V> dahu_distance_color(tree_t T, image2d<V> ima_compo )
    {
        // Initiation

        unsigned nodenumbers = T.nodes().size();

        auto& nodes1 = T._get_data()->m_nodes;
        auto& S     = T._get_data()->m_S;
        auto& pmap  = T._get_data()->m_pmap;

        box2d D = pmap.domain();

        image2d <V>  distance_map(D);
        extension::fill(distance_map, 0);

        image2d <uint8_t>  distance_map1(D);
        extension::fill(distance_map1, 0);

        std::vector<V>  min_node(nodenumbers);
        std::vector<V>  max_node(nodenumbers);

        std::vector<V>  min_temp_node(nodenumbers);
        std::vector<V>  max_temp_node(nodenumbers);

        std::vector<V>  dist(nodenumbers);
        std::vector<unsigned>  dist_gray_temp_node(nodenumbers);


        std::vector< std::vector<unsigned> > children(nodenumbers);

        std::vector<int>  seed_node(nodenumbers);
        std::fill (seed_node.begin(),seed_node.end(),0);
        seed_node[0] = 1;
        std::vector<unsigned>  dejavu(nodenumbers);
        std::fill (dejavu.begin(),dejavu.end(),0);

        unsigned height = ima_compo.nrows();
        unsigned width = ima_compo.ncols();


        mln_foreach(auto p, D)
        {
            if (is_border(p,height,width))
            {
                dejavu[pmap(p)] = 1;
                seed_node[pmap(p)] = 1;
            }
        }






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

                min_temp_node[x.id()][i] = medianvalue1[i];
                max_temp_node[x.id()][i] = medianvalue1[i];
                dist_gray_temp_node[x.id()] = 255 *3;


            }


            unsigned q = x.get_parent_id();
            children[q].push_back(x.id());

        }

        // chay up len


        std::cout << "chay up len  "  << std::endl;

        mln_reverse_foreach(auto x, T.nodes())
        {
            unsigned dis_temp = 0;
            V min_temp = min_node[x.get_parent_id()];
            V max_temp = max_node[x.get_parent_id()];

            if (dejavu[x.id()] == 1)
            {
                unsigned par_p = x.get_parent_id();

                if (seed_node[par_p] == 0)
                {

                    for (int i = 0; i< 3; ++i)
                    {
                        if (min_node[par_p][i] > min_temp_node[x.id()][i])
                            min_temp[i] = min_temp_node[x.id()][i];

                        if (max_node[par_p][i] < max_temp_node[x.id()][i])
                            max_temp[i] = max_temp_node[x.id()][i];

                        dis_temp = dis_temp + max_temp[i]  - min_temp[i];

                    }

                    if (dis_temp < dist_gray_temp_node[par_p])
                    {
                        dist_gray_temp_node[par_p] = dis_temp;
                        min_temp_node[par_p] = min_temp;
                        max_temp_node[par_p] = max_temp;
                    }
                    dejavu[par_p] =1;

                }
            }
        }

        // chay down xuong


        std::cout << "chay down xuong "  << std::endl;

        std::cout << "root " << T.get_root_id()  << std::endl;

        mln_foreach (auto x, T.nodes())
        {
            V min_temp = min_node[x.id()];
            V max_temp = max_node[x.id()];
            unsigned dis_temp = 0;

            if(seed_node[x.id()] == 0)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if (min_node[x.id()][i] > min_temp_node[x.get_parent_id()][i])
                        min_temp[i] = min_temp_node[x.get_parent_id()][i];
                    if (max_node[x.id()][i] < max_temp_node[x.get_parent_id()][i])
                        max_temp[i] = max_temp_node[x.get_parent_id()][i];
                    dis_temp = dis_temp + max_temp[i]  - min_temp[i];
                }

                if (dis_temp < dist_gray_temp_node[x.id()])
                {
                    dist_gray_temp_node[x.id()] = dis_temp;
                    max_temp_node[x.id()] = max_temp;
                    min_temp_node[x.id()] = min_temp;
                }
            }
        }


        // back propagation

        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                distance_map(ima_compo.point_at_index(p)) = max_temp_node[x.id()] - min_temp_node[x.id()];
            }
        }
        return distance_map;

    }


    image2d<uint8> dahu_distance(tree_t T, image2d<float> U )

    {

        // Initiation

        unsigned nodenumbers = T.nodes().size();

        auto& nodes1 = T._get_data()->m_nodes;
        auto& S     = T._get_data()->m_S;
        auto& pmap  = T._get_data()->m_pmap;

        box2d D = pmap.domain();
        std::cout <<"Domain  " <<D << std::endl;
        unsigned height = D.pmax[0];
        unsigned width = D.pmax[1];

        image2d <uint8>  distance_map(D);
        extension::fill(distance_map, 0);

        std::vector<float>  min_node(nodenumbers);
        std::vector<float>  max_node(nodenumbers);

        std::vector<float>  min_temp_node(nodenumbers);
        std::vector<float>  max_temp_node(nodenumbers);

        //std::vector<V>  dist(nodenumbers);
        std::vector<unsigned>  dist_gray_temp_node(nodenumbers);



        std::vector<int>  seed_node(nodenumbers);
        std::fill (seed_node.begin(),seed_node.end(),0);
        seed_node[0] = 1;
        std::vector<unsigned>  dejavu(nodenumbers);
        std::fill (dejavu.begin(),dejavu.end(),0);


        point2d seedpoints;
        seedpoints = point2d (69,372);

        std::cout << "seed "  << pmap(seedpoints)  << std::endl;


        mln_foreach(auto p, D)
        {
            if (is_border(p,height,width))
            {
                dejavu[pmap(p)] = 1;
                seed_node[pmap(p)] = 1;
            }
        }




        mln_foreach(auto x, T.nodes())
        {

            unsigned p = x.first_point();


            min_node[x.id()] = U(U.point_at_index(p));
            min_temp_node[x.id()] = U(U.point_at_index(p));
            max_node[x.id()] = U(U.point_at_index(p));
            max_temp_node[x.id()] = U(U.point_at_index(p));
            dist_gray_temp_node[x.id()] = INT_MAX;


        }

        // chay up len


        std::cout << "chay up len  "  << std::endl;

        mln_reverse_foreach(auto x, T.nodes())
        {
            unsigned dis_temp = 0;
            float min_temp = min_node[x.get_parent_id()];
            float max_temp = max_node[x.get_parent_id()];

            if (dejavu[x.id()] == 1)
            {
                unsigned par_p = x.get_parent_id();

                if (seed_node[par_p] == 0)
                {

                    if (min_node[par_p] > min_temp_node[x.id()])
                        min_temp = min_temp_node[x.id()];

                    if (max_node[par_p] < max_temp_node[x.id()])
                        max_temp = max_temp_node[x.id()];

                    dis_temp = dis_temp + max_temp  - min_temp;


                    if (dis_temp < dist_gray_temp_node[par_p])
                    {
                        dist_gray_temp_node[par_p] = dis_temp;
                        min_temp_node[par_p] = min_temp;
                        max_temp_node[par_p] = max_temp;
                    }
                    dejavu[par_p] =1;

                }
            }
        }

        // chay down xuong


        std::cout << "chay down xuong "  << std::endl;


        mln_foreach (auto x, T.nodes())
        {
            float min_temp = min_node[x.id()];
            float max_temp = max_node[x.id()];
            unsigned dis_temp = 0;

            if(seed_node[x.id()] == 0)
            {

                if (min_node[x.id()] > min_temp_node[x.get_parent_id()])
                    min_temp = min_temp_node[x.get_parent_id()];
                if (max_node[x.id()] < max_temp_node[x.get_parent_id()])
                    max_temp = max_temp_node[x.get_parent_id()];
                dis_temp = dis_temp + max_temp  - min_temp;


                if (dis_temp < dist_gray_temp_node[x.id()])
                {
                    dist_gray_temp_node[x.id()] = dis_temp;
                    max_temp_node[x.id()] = max_temp;
                    min_temp_node[x.id()] = min_temp;
                }
            }
        }


        // back propagation

        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                distance_map(U.point_at_index(p)) = max_temp_node[x.id()] - min_temp_node[x.id()];
            }
        }
        return distance_map;
    }

}



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] α₀ α₁ λ output[rgb]  output[rgb]\n"
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
  int a0 = std::atoi(argv[2]);
  int a1 = std::atoi(argv[3]);
  int lambda = std::atoi(argv[4]);
  const char* output_path_MTOS = argv[5];
  const char* output_path_TOS = argv[6];


  tbb::task_scheduler_init init;


  std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_resized";
  std::string outputDirectory_mtos = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/MTOS";
  std::string outputDirectory_mtos_scalar = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/MTOS_scalar";
  std::string outputDirectory_fusion = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/fusion";
  std::string outputDirectory_r = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/R";
  std::string outputDirectory_g = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/G";
  std::string outputDirectory_b = "/home/minh/Documents/code/code_edwin/build/apps/papier/MTOS_vs_3chanels/B";

  //std::cout << inputDirectory << std::endl;
  DIR *directory = opendir (inputDirectory.c_str());
//  DIR *directory1 = opendir (outputDirectory.c_str());
  struct dirent *_dirent = NULL;
  if(directory == NULL)
  {
      printf("Cannot open Input Folder\n");
      return 1;
  }

//  if(directory1 == NULL)
//  {
//      printf("Cannot open Output Folder\n");
//      return 1;
//  }

  while((_dirent = readdir(directory)) != NULL)
  {

    if ((std::string(_dirent->d_name) != ".") and (std::string(_dirent->d_name) != "..") )
    {


          std::string fileName = inputDirectory + "/" +std::string(_dirent->d_name);

          // 1. Compute the individual ToS
          using namespace mln;
          typedef rgb8 V;

          image2d<V> ima;
          io::imread(fileName.c_str(), ima);

          image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
          image2d<V> F = interpolate_k1(f);


          //std::string outputDirectory = "/lrde/home/movn/Documents/code/code_edwin2/program/apps/mumford_shah_on_tree/abc";

          image2d<V> ima_compo(F.domain() );



          // COmpute the color tree of shapes T
          tree_t T;

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

            t[0]._reorder_pset();
            t[1]._reorder_pset();
            t[2]._reorder_pset();




            std::ofstream file("abc.dot",
                     std::ios_base::out | std::ios_base::trunc);
            write_graphviz(file, T);





            image2d<rgb8 > dahu_MTOS = dahu_distance_color( T,  ima_compo );
            image2d<uint8> dahu_MTOS_scalar(F.domain());

            image2d<uint8 > dahu_TOS_0 = dahu_distance( t[0],  U0 );
            image2d<uint8 > dahu_TOS_1 = dahu_distance( t[1],  U1 );
            image2d<uint8 > dahu_TOS_2 = dahu_distance( t[2],  U2 );

            unsigned alpha = 20;
            unsigned beta = 210;
            int max1 = 0;
            int maxr = 0;
            int maxg = 0;
            int maxb = 0;



            mln_foreach(auto p, F.domain())
            {
                dahu_MTOS_scalar(p) = dahu_MTOS(p)[0]/3 + dahu_MTOS(p)[1]/3 + dahu_MTOS(p)[2]/3;
                if (dahu_MTOS_scalar(p) > max1)
                    max1 = dahu_MTOS_scalar(p);
                if (dahu_TOS_0(p) > maxr)
                    maxr = dahu_TOS_0(p);
                if (dahu_TOS_1(p) > maxg)
                    maxg = dahu_TOS_1(p);
                if (dahu_TOS_2(p) > maxb)
                    maxb = dahu_TOS_2(p);

            }

            mln_foreach(auto p , F.domain())
            {
                dahu_MTOS_scalar(p) = uint8(float(dahu_MTOS_scalar(p))/ float(max1) *255);
                dahu_TOS_0(p) = uint8(float(dahu_TOS_0(p))/ float(maxr) *255);
                dahu_TOS_1(p) = uint8(float(dahu_TOS_1(p))/ float(maxg) *255);
                dahu_TOS_2(p) = uint8(float(dahu_TOS_2(p))/ float(maxb) *255);

            }

            mln_foreach(auto p , F.domain())
            {
                if (dahu_MTOS_scalar((p)) <= alpha)
                    dahu_MTOS_scalar((p)) = 0;
                else if (dahu_MTOS_scalar((p)) >= beta)
                    dahu_MTOS_scalar((p)) = 255;
                else
                {
                    dahu_MTOS_scalar((p)) = uint8(float(dahu_MTOS_scalar((p)) - alpha)/ float(beta - alpha) *255);
                }

                if (dahu_TOS_0((p)) <= alpha)
                    dahu_TOS_0((p)) = 0;
                else if (dahu_TOS_0((p)) >= beta)
                    dahu_TOS_0((p)) = 255;
                else
                {
                    dahu_TOS_0((p)) = uint8(float(dahu_TOS_0((p)) - alpha)/ float(beta - alpha) *255);
                }

                if (dahu_TOS_1((p)) <= alpha)
                    dahu_TOS_1((p)) = 0;
                else if (dahu_TOS_1((p)) >= beta)
                    dahu_TOS_1((p)) = 255;
                else
                {
                    dahu_TOS_1((p)) = uint8(float(dahu_TOS_1((p)) - alpha)/ float(beta - alpha) *255);
                }

                if (dahu_TOS_2((p)) <= alpha)
                    dahu_TOS_2((p)) = 0;
                else if (dahu_TOS_2((p)) >= beta)
                    dahu_TOS_2((p)) = 255;
                else
                {
                    dahu_TOS_2((p)) = uint8(float(dahu_TOS_2((p)) - alpha)/ float(beta - alpha) *255);
                }

            }





            fileName = outputDirectory_mtos + "/" + std::string(_dirent->d_name);
            io::imsave(dahu_MTOS, fileName.c_str());

            fileName = outputDirectory_mtos_scalar + "/" + std::string(_dirent->d_name);
            io::imsave(dahu_MTOS_scalar, fileName.c_str());

//            fileName = outputDirectory_fusion + "/" + std::string(_dirent->d_name);
//            io::imsave(dahu_TOS, fileName.c_str());

            fileName = outputDirectory_r + "/" + std::string(_dirent->d_name);
            io::imsave(dahu_TOS_0, fileName.c_str());

            fileName = outputDirectory_g + "/" + std::string(_dirent->d_name);
            io::imsave(dahu_TOS_1, fileName.c_str());

            fileName = outputDirectory_b + "/" + std::string(_dirent->d_name);
            io::imsave(dahu_TOS_2, fileName.c_str());

        }
    }

    closedir(directory);

}

