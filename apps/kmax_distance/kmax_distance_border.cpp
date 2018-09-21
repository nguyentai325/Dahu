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
#include <numeric>



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
        if (p[0] == 0 or p[0]== height -1 or p[1] == 0 or p[1] == width)
            return true;
        else
            return false;

    }

    bool is_border1 (point2d p, unsigned height , unsigned width)
    {
        if (p[0] != 0  and p[0]  != 1  and p[0] != height -1 and p[0] != height -2 and p[1] != 0  and p[1]  != 1 and p[1]  != width - 1 and p[1] != width -2)
        {
            if (p[0] == 2 or p[0] == height - 3 or p[1] == width -3 or p[1] == 2  )
                return true;
            else
                return false;
        }
    }

    image2d<uint8> rgb2gray(const image2d<rgb8>& input)
    {
        image2d<uint8> output(input.nrows() ,
                            input.ncols() );
        box2d dom = output.domain();


        mln_foreach(auto p,dom)
        {
            output(p) = 0.2989 * input(p)[0] + 0.5870 * input(p)[1] + 0.1140 * input(p)[2];
        }

        return output;
    }

    image2d<uint8> addborder_gray(const image2d<uint8> ima)
    {
        image2d<uint8> out(ima.nrows() + 2, ima.ncols() + 2);

        {
          box2d box = ima.domain();
          box.pmin += 1; box.pmax += 1;
          copy(ima, out | box);
        }

        uint8 median;
        unsigned ncols = ima.ncols(), nrows = ima.nrows();
        {
          std::vector<uint8> border;
          border.reserve(2 * (nrows + ncols) - 4);
          for (unsigned i = 0; i < ncols; ++i) {
        border.push_back(ima.at(0,i));
        border.push_back(ima.at(nrows-1,i));
          }

          for (unsigned i = 1; i < nrows-1; ++i) {
        border.push_back(ima.at(i,0));
        border.push_back(ima.at(i,ncols-1));
          }

          std::sort(border.begin(),  border.end(), std::less<uint8>());
          if (border.size() % 2 == 0) {
        //V a = border[border.size()/2 - 1], b = border[border.size()/2];
        //median = a + (b-a) / 2;
        median = border[border.size()/2];
          } else
        median = border[border.size()+1/2];
        }

        {
          for (unsigned i = 0; i < ncols+2; ++i) {
        out.at(0,i) = median;
        out.at(nrows+1,i) = median;
          }

          for (unsigned i = 1; i < nrows+1; ++i) {
        out.at(i,0) = median;
        out.at(i,ncols+1) = median;
          }
        }
        return out;
    }

    image2d<uint8> interpolate(const image2d<uint8>& ima)
    {
      image2d<uint8> out(2*ima.nrows()-1, 2*ima.ncols()-1);
      typedef point2d P;
      mln_foreach(point2d p, ima.domain())
        {
      uint8 a = ima.at(p),
        b = ima.at(p + P{0,1}),
        c = ima.at(p + P{1,0}),
        d = ima.at(p + P{1,1});

      point2d q = 2 * p;
      out.at(q) = ima.at(p);
      out.at(q + P{0,1}) = (a + b) / 2;
      out.at(q + P{1,0}) = (a + c) / 2;
      out.at(q + P{1,1}) = (a + b + c + d) / 4;
        }

      return out;
    }
}



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] α₀ α₁ λ output[rgb]\n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}


int main(int argc, char** argv)
{
  if (argc < 5)
    usage(argv);



  const char* input_path = argv[1];
  int a0 = std::atoi(argv[2]);
  int a1 = std::atoi(argv[3]);
  int lambda = std::atoi(argv[4]);
  const char* output_path = argv[5];


  tbb::task_scheduler_init init;


  std::string inputDirectory = "/media/minh/DATA/Study/dataset/test_resized";
  std::string outputDirectory = "/media/minh/DATA/Study/Results/K_max/K_max100";
  //std::cout << inputDirectory << std::endl;
  DIR *directory = opendir (inputDirectory.c_str());
  DIR *directory1 = opendir (outputDirectory.c_str());
  struct dirent *_dirent = NULL;
  if(directory == NULL)
  {
      printf("Cannot open Input Folder\n");
      return 1;
  }

  if(directory1 == NULL)
  {
      printf("Cannot open Output Folder\n");
      return 1;
  }


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

          image2d<uint8>  ima_gray = rgb2gray(ima);






          image2d<uint8> f = addborder_gray(ima_gray);

          image2d<uint8> F = interpolate(f);


          // build a tree of shape of gray image

          tree_t T;
          T = morpho::cToS(f, c4);
          if (a0 > 0)
          {
            grain_filter_inplace(T, a0);
            T.shrink_to_fit();
          }
          T._reorder_pset();



          auto& U  = T._get_data()->m_Uv;

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

          std::vector<int>  seed_node(nodenumbers);
          std::fill (seed_node.begin(),seed_node.end(),0);
          seed_node[0] = 1;
          std::vector<unsigned>  dejavu(nodenumbers);
          std::fill (dejavu.begin(),dejavu.end(),0);

          // put seeds on the boundary of image

          std::vector<uint8>  value_node(nodenumbers);
          //const unsigned k = 2;
          typedef std::array<uint8, 200> fixsize;
          std::vector<fixsize> k_max(nodenumbers);


          mln_foreach(auto p, D)
          {
              if (is_border(p,height,width))
              {
                  dejavu[pmap(p)] = 1;
                  seed_node[pmap(p)] = 1;

              }
          }




          //// Node value for each node on the tree

//          mln_foreach(auto x, T.nodes())
//          {
//              unsigned p = x.first_point();
//              value_node[x.id()] = U(U.point_at_index(p));
//              k_max[x.id()].fill(0);
//          }

          mln_foreach(auto x, T.nodes())
          {

              std::vector<uint8> listpixel1;

              mln_foreach (auto p, x.proper_pset())
              {
                  listpixel1.push_back(U(U.point_at_index(p)));
              }

              std::sort(listpixel1.begin(), listpixel1.end());
              value_node[x.id()] = listpixel1[listpixel1.size()/2];
              k_max[x.id()].fill(0);

          }


          // chay up len


          std::cout << "chay up len  "  << std::endl;

          mln_reverse_foreach(auto x, T.nodes())
          {


              if (dejavu[x.id()] == 1)
              {
                  unsigned par_p = x.get_parent_id();
                  //std::sort(k_max[x.id()].begin(), k_max[x.id()].end(), std::greater<uint8>());
                  k_max[par_p] = k_max[x.id()];

                  if (seed_node[par_p] == 0)
                  {
                      uint8 distance_temp = std::abs(value_node[par_p]- value_node[x.id()]) ;

                      if (distance_temp > k_max[par_p][199])
                          k_max[par_p][199] = distance_temp;
                      std::sort(k_max[par_p].begin(), k_max[par_p].end(), std::greater<uint8>());

                      dejavu[par_p] =1;

                  }
              }
          }

          // chay down xuong


          std::cout << "chay down xuong "  << std::endl;
          auto root_id = T.get_root();


          mln_foreach (auto x, T.nodes())
          {


              if(seed_node[x.id()] == 0 and x.id() != T.get_root_id())
              {


                  k_max[x.id()] = k_max[x.get_parent_id()];




                  uint8 distance_temp = std::abs(value_node[x.get_parent_id()]- value_node[x.id()]) ;
                  //std::cout << "distance "  << int(distance_temp) << std::endl;
                  if (distance_temp > k_max[x.id()][199])
                      k_max[x.id()][199] = distance_temp;

//                  for (unsigned i = 0 ; i < 10 ; i++)
//                        std::cout  <<int(k_max[x.id()][i])  << std::endl;

                  std::sort(k_max[x.id()].begin(), k_max[x.id()].end(), std::greater<uint8>());

//                  for (unsigned i = 0 ; i < 10 ; i++)
//                        std::cout << "sort  " <<int(k_max[x.id()][i])  << std::endl;
              }
          }


          point2d pp = point2d(height/2,width/2);

          // back propagation

          unsigned alpha = 20;
          unsigned beta = 160;

          mln_foreach(auto x, T.nodes())
          {
              mln_foreach (auto p, x.proper_pset())
              {
                  uint8 sum_of_elems = 0;
//                  std::for_each(k_max[x.id()].begin(), k_max[x.id()].end(), [&] (uint8 n) {
//                      sum_of_elems += n;});
                  for (int t = 0; t < 200 ; t++)
                  {
                      sum_of_elems = sum_of_elems + k_max[x.id()][t];
                  }
                  distance_map(F.point_at_index(p)) =  sum_of_elems;
                  //std::cout << int(distance_map(F.point_at_index(p)))  << std::endl;

//                  if (distance_map(F.point_at_index(p)) <= alpha)
//                      distance_map(F.point_at_index(p)) = 0;
//                  else if (distance_map(F.point_at_index(p)) >= beta)
//                      distance_map(F.point_at_index(p)) = 255;
//                  else
//                  {
//                      distance_map(F.point_at_index(p)) = uint8(float(distance_map(F.point_at_index(p)) - alpha)/ float(beta - alpha) *255);
//                  }
              }
          }



          fileName = outputDirectory + "/" + std::string(_dirent->d_name);

          io::imsave(distance_map, fileName.c_str());

      }
  }

  closedir(directory);



}

