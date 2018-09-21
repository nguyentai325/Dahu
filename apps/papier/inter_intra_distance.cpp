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
#include <mln/morpho/structural/dilate.hpp>

#include "immerse.hpp"

#include "function.hpp"





#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"

#include <chrono>


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

  using namespace mln;
  typedef rgb8 V;

  const char* input_path = argv[1];
  const char* input_path_gt = argv[2];
  int a0 = std::atoi(argv[3]);
  int a1 = std::atoi(argv[4]);
  int lambda = std::atoi(argv[5]);
  const char* output_path = argv[6];


  tbb::task_scheduler_init init;


  std::string inputDirectory = "/media/minh/DATA/Study/database/MSRA10K_Imgs_GT/MSRA10K_Imgs_GT/Imgs";
  std::string inputDirectory1 = "/media/minh/DATA/Study/database/MSRA10K_Imgs_GT/MSRA10K_Imgs_GT/Imgs";


  //std::cout << inputDirectory << std::endl;
  DIR *directory = opendir (inputDirectory.c_str());
  struct dirent *_dirent = NULL;
  if(directory == NULL)
  {
      printf("Cannot open Input Folder\n");
      return 1;
  }


  std::vector<uint8> inter_distance_dahu;
  double inter_distance_dahu_sum = 0;
  std::vector<uint8> intra_distance_dahu;
  double intra_distance_dahu_sum = 0;


  std::vector<uint8> inter_distance_water;
  std::vector<uint8> intra_distance_water;
  double inter_distance_water_sum = 0;
  double intra_distance_water_sum = 0;


  std::vector<uint8> inter_distance_mst;
  std::vector<uint8> intra_distance_mst;
  double inter_distance_mst_sum = 0;
  double intra_distance_mst_sum = 0;


  while((_dirent = readdir(directory)) != NULL)
  {

    if ((std::string(_dirent->d_name) != ".") and (std::string(_dirent->d_name) != "..") )
    {


          std::string fileName = inputDirectory + "/" +std::string(_dirent->d_name);
          std::string name = std::string(_dirent->d_name);


          std::istringstream iss(name);
          std::vector<std::string> tokens;
          std::string token;
          while (std::getline(iss, token, '.')) {
              if (!token.empty())
                  tokens.push_back(token);
          }


          std::string surname = "jpg";

          if (surname.compare(tokens[1]) == 0)
          {

          //std::string fileName = input_path;
          std::string fileName_gt = inputDirectory1 + "/" +tokens[0] +".png";
          //std::cout << token << std::endl;


          //std::string fileName = input_path;

          // 1. Compute the individual ToS
          std::cout << fileName  << std::endl;
          std::cout << fileName_gt  << std::endl;


          image2d<V> ima;
          io::imread(fileName.c_str(), ima);

          image2d<uint8> ima_gt;
          io::imread(fileName_gt.c_str(), ima_gt);

            //image2d<uint8>  ima_gt = rgb2gray(ima_gt1);
          //////////////////////////////////////  Dahu distance /////////////////////////////////////////////////////////////////

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


          // Initiation

          unsigned nodenumbers = T.nodes().size();

          auto& nodes1 = T._get_data()->m_nodes;
          auto& S     = T._get_data()->m_S;
          auto& pmap  = T._get_data()->m_pmap;

          box2d D = pmap.domain();

          std::vector<V>  node_value(nodenumbers);


          // node value

          mln_foreach(auto x, T.nodes())
          {
              std::vector<V> listpixel1;
              mln_foreach (auto p, x.proper_pset())
              {
                  listpixel1.push_back(ima_compo(ima_compo.point_at_index(p)));
              }
              std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end(),lexicographicalorder_less<value_t>());
              V medianvalue1;
              node_value[x.id()] = listpixel1[listpixel1.size()/2];
          }

          // depth image

          std::vector<unsigned>  depth_node = depth_node_compute(T);


          const unsigned N = 20;
          box2d Dom = f.domain();

          std::vector<point2d> pt(N);
          std::vector<bool> label(N);










          for (unsigned i = 0; i < N; ++i)
          {
            pt[i] = point2d(std::rand() % Dom.pmax[0], std::rand() % Dom.pmax[1]);
            if (ima_gt(pt[i]) == 0)
                label[i] = false;
            else
                label[i] = true;
          }






          auto t1 = std::chrono::high_resolution_clock::now();
          uint8 d_dahu;
          for (unsigned i = 1; i < N; i += 1)
          {
                d_dahu = compute_dahu(T, depth_node, node_value,  pt[0], pt[i]);
                if (label[i] == label[0])
                {
                    intra_distance_dahu.push_back(d_dahu);
                    intra_distance_dahu_sum = intra_distance_dahu_sum + d_dahu;
                }
                else
                {
                    inter_distance_dahu.push_back(d_dahu);
                    inter_distance_dahu_sum = inter_distance_dahu_sum + d_dahu;
                }
          }

          //std::cout << int(d_dahu)  << std::endl;

          auto t2 = std::chrono::high_resolution_clock::now();

          std::chrono::duration<double,std::milli> elapsed = t2 - t1;

          //std::cout << "Algorithm Runtime is: " << elapsed.count() << " miliseconds." << std::endl;

        //  image2d<rgb8>  dahu_image = compute_dahu_image(T,ima_compo);
        //  io::imsave(dahu_image, "dmap_dahu.png");






          //////////////////////////////////////  Waterflow MBD /////////////////////////////////////////////////////////////////






          t1 = std::chrono::high_resolution_clock::now();
          uint8 d_waterflow;
          for (unsigned i = 1; i < N; i += 1)
          {
            d_waterflow = compute_waterflow( f, pt[0], pt[i]);
            if (label[i] == label[0])
            {
                intra_distance_water.push_back(d_waterflow);
                intra_distance_water_sum = intra_distance_water_sum + d_waterflow;
            }
            else
            {
                inter_distance_water.push_back(d_waterflow);
                inter_distance_water_sum = inter_distance_water_sum + d_waterflow;
            }
          }

          //std::cout << int(d_waterflow)   << std::endl;

          t2 = std::chrono::high_resolution_clock::now();

          elapsed = t2 - t1;

          //std::cout << "Algorithm Runtime is: " << elapsed.count() << " miliseconds." << std::endl;


        //  image2d<rgb8>  waterflow_image = waterflow_mbd_distance(f);
        //  io::imsave(waterflow_image, "dmap_waterflow.png");





          //////////////////////////////////////  MST MBD /////////////////////////////////////////////////////////////////


          std::vector<point2d> S1;
          image2d<point2d> parent = primMST_color(f, S1, point2d(0,0));


          image2d<unsigned> depth_node_mst = depth_node_mst_compute(parent, S1 );





          t1 = std::chrono::high_resolution_clock::now();
          uint8 d_mstmbd;
          for (unsigned i = 1; i < N; i += 1)
          {
            d_mstmbd = compute_mstmbd( f, parent, depth_node_mst, pt[0], pt[i]);
            if (label[i] == label[0])
            {
                intra_distance_mst.push_back(d_mstmbd);
                intra_distance_mst_sum = intra_distance_mst_sum + d_mstmbd;
            }
            else
            {
                inter_distance_mst.push_back(d_mstmbd);
                inter_distance_mst_sum = inter_distance_mst_sum + d_mstmbd;
            }
          }

          //std::cout << int(d_mstmbd)   << std::endl;

          t2 = std::chrono::high_resolution_clock::now();

          elapsed = t2 - t1;

          //std::cout << "Algorithm Runtime is: " << elapsed.count() << " miliseconds." << std::endl;


        //  image2d<rgb8>  waterflow_image = waterflow_mbd_distance(f);
        //  io::imsave(waterflow_image, "dmap_waterflow.png");

        //  std::cout << inter_distance_mst_sum  << "  "  << intra_distance_mst_sum  << std::endl;
        //  std::cout << intra_distance_mst.size()  << "  "  << inter_distance_mst.size()  << std::endl;
          }




        }
    }

    closedir(directory);

    std::cout << (float(inter_distance_water_sum)*intra_distance_water.size())/(inter_distance_water.size()*float(intra_distance_water_sum))  << "  "  << (float(inter_distance_dahu_sum)*intra_distance_dahu.size())/(float(inter_distance_dahu.size()*intra_distance_dahu_sum))  << std::endl;

    std::cout << (float(inter_distance_mst_sum)*intra_distance_mst.size())/(inter_distance_mst.size()*float(intra_distance_mst_sum))  << "  "  << (float(inter_distance_dahu_sum)*intra_distance_dahu.size())/(float(inter_distance_dahu.size()*intra_distance_dahu_sum))  << std::endl;



}

