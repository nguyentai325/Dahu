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
#include<bits/stdc++.h>

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



    // create a minimum spanning tree with Kruskal algorithm
    typedef  std::pair<int, int> iPair;

    struct GraphMST
    {
        int V, E;
        std::vector< std::pair<int, iPair> > edges;
        std::vector<int> zpar;
        std::vector<int> parent;
        std::vector<int> rnk;
        std::vector<int> level;
        std::vector<int> nodes;
        //std::vector< std::vector<int> >  children;
        //children.reserve(17);


        // Constructor
        GraphMST(int V, int E)
        {
            this->V = V;
            this->E = E;
            //children.reserve(2*V-1);
        }

        // Utility function to add an edge

        void initiate(int V)
        {
            for (int i = 0; i < V; i++)
            {
                rnk.push_back(0);
                parent.push_back(i);
                zpar.push_back(i);
                level.push_back(0);
                nodes.push_back(i);
                //children.push_back(std::vector<int>());
            }

        }
        void addEdge(int u, int v, int w)
        {
            edges.push_back({w, {u, v}});
        }

        // Function to find MST using Kruskal's
        // MST algorithm


        int find(int u)
        {
            /* Make the parent of the nodes in the path
               from u--> parent[u] point to parent[u] */
            if (u != zpar[u])
                zpar[u] = find(zpar[u]);
            return zpar[u];
        }

        // Union by rank
        void merge(int x, int y, int dt, int numbernodes)
        {
            //x = find(x), y = find(y);


            /* Make tree with smaller height
               a subtree of the other tree  */
            nodes.push_back(numbernodes);
            parent.push_back(numbernodes);
            zpar.push_back(numbernodes);

            rnk.push_back(0);
            parent[x] = numbernodes;
            parent[y] = numbernodes;
            zpar[x] = numbernodes;
            zpar[y] = numbernodes;
            //children[numbernodes].push_back(x);
            //children[numbernodes].push_back(y);
            std::cout << x  <<  "    "   << y << "   has parent    "   << numbernodes  << std::endl;
            rnk[x] = rnk[x] + 1;
            rnk[y] = rnk[y] + 1;
            level.push_back(dt);

        }

        int kruskalMST();

        void compute_rnk ();

        int LCA(int x, int y);

        image2d<uint8> saliencymap();


    };



     /* Functions returns weight of the MST*/

    int GraphMST::kruskalMST()
    {
        int mst_wt = 0; // Initialize result

        // Sort edges in increasing order on basis of cost
        std::sort(edges.begin(), edges.end());

        // Create disjoint sets
        initiate(V);

        // Iterate through all sorted edges
        int numbernodes = V;
        std::vector< std::pair<int, iPair> >::iterator it;
        for (it=edges.begin(); it!=edges.end(); it++)
        {
            int u = it->second.first;
            int v = it->second.second;

            int set_u = find(u);
            int set_v = find(v);

            int dt = it->first;

            // Check if the selected edge is creating
            // a cycle or not (Cycle is created if u
            // and v belong to same set)
            if (set_u != set_v)
            {

                // Current edge will be in the MST

                // Update MST weight
                mst_wt += it->first;

                // Merge two sets
                merge(set_u, set_v,dt,numbernodes);
                numbernodes = numbernodes + 1;
            }
        }

        return mst_wt;
    }

    void GraphMST::compute_rnk()
    {

        for (int i = nodes.size() -2  ; i >= 0  ; i--)
        {
            rnk[nodes[i]] = rnk[parent[nodes[i]]] + 1;
            std::cout << i <<   "   "  << nodes[i]  << "   "  << rnk[nodes[i]] << std::endl;
        }

    }


    int GraphMST::LCA(int x, int y)
    {

        int node_cur = rnk[x] >= rnk[y] ? x : y;
        int node_obj = rnk[x] >= rnk[y] ? y :x  ;

        int d_cur = rnk[node_cur];
        int d_obj = rnk[node_obj];


        while (d_cur > d_obj)
        {
            d_cur -= 1;
            //path.push_back(node_cur);
            node_cur = parent[node_cur];
        }
        while (node_cur != node_obj)
        {
            //path.push_back(node_cur);
            //path.push_back(node_obj);
            node_cur = parent[node_cur];
            node_obj = parent[node_obj];
        }

        //path.push_back(node_cur);

        int lca = node_cur;
        //cout << "lca   " <<  lca  << endl;

        return lca;

    }



    image2d<uint8> GraphMST::saliencymap()
    {

        image2d<uint8> saliency(17, 17);
        box2d D1 = saliency.domain();
        mln_foreach(auto p , D1)
        {
            saliency(p) = 0;
        }

        std::cout << "saliency map "  << std::endl;
        std::vector< std::pair<int, iPair> >::iterator it;
        for (it=edges.begin(); it!=edges.end(); it++)
        {
            int u = it->second.first;
            int v = it->second.second;
            int lca1 = LCA(u,v);
            std::cout <<  u  << "   "<< v  << "  "<< lca1  << "   "<<   level[lca1]  << std::endl;
            point2d point_u = point2d(u/9,u%9);
            point2d point_v = point2d(v/9,v%9);
            saliency(point_u+point_v) =level[lca1];
        }

        mln_foreach(auto p , D1)
        {
            if (p[0]%2 == 1 and p[1]%2 ==1)
            {
                saliency(p) = std::max(saliency(point2d(p[0]-1,p[1])),std::max(saliency(point2d(p[0]+1,p[1])),std::max(saliency(point2d(p[0],p[1]-1)),saliency(point2d(p[0],p[1]+1)))));
            }
        }

        return saliency;


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


//  while((_dirent = readdir(directory)) != NULL)
//  {

//    if ((std::string(_dirent->d_name) != ".") and (std::string(_dirent->d_name) != "..") )
//    {


//          std::string fileName = inputDirectory + "/" +std::string(_dirent->d_name);
//          std::string name = std::string(_dirent->d_name);


//          std::istringstream iss(name);
//          std::vector<std::string> tokens;
//          std::string token;
//          while (std::getline(iss, token, '.')) {
//              if (!token.empty())
//                  tokens.push_back(token);
//          }


//          std::string surname = "jpg";

//          if (surname.compare(tokens[1]) == 0)
//          {

//              //std::string fileName = input_path;
//              std::string fileName_gt = inputDirectory1 + "/" +tokens[0] +".png";
//              //std::cout << token << std::endl;


//              //std::string fileName = input_path;

//              // 1. Compute the individual ToS
//              std::cout << fileName  << std::endl;
//              std::cout << fileName_gt  << std::endl;


//              image2d<V> ima;
//              io::imread(fileName.c_str(), ima);

//              image2d<uint8> ima_gt;
//              io::imread(fileName_gt.c_str(), ima_gt);

                //image2d<uint8>  ima_gt = rgb2gray(ima_gt1);
              //////////////////////////////////////  Dahu distance /////////////////////////////////////////////////////////////////


              std::string fileName = input_path;
              image2d<V> ima;
              io::imread(fileName.c_str(),ima);
              image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
              image2d<V> F = interpolate_k1(f);
              unsigned height = f.nrows();
              unsigned width = f.ncols();



              image2d<V> image(f.domain());

              mln_foreach( auto p, f.domain())
              {
                  image(p) = f(p);
              }

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

              unsigned Nx = 10;
              unsigned Ny = 10;


              // create a set of points

              point2d pointxy[9][9];


              for (int i= 0 ; i < Ny-1 ; i++)
              {
                  int pointy = height/Ny*(i+1);
                  for (int j = 0; j < Nx-1 ; j++)
                  {
                      int pointx = width/Nx*(j+1);
                      pointxy[i][j] = point2d(pointy,pointx);
                      image(pointxy[i][j]) = {255,0,0};
                      //std::cout << pointxy[i][j]  << std::endl;

                  }
              }

              io::imsave(image, "dm.png");


              // compute a distance between a set of point to construct a minimum spanning tree
              // create a MST of these points


              int Vertex =81, Edge=9*8*2;
              GraphMST g(Vertex, Edge);




              int dy[2] = {0, 1};
              int dx[2] = {1, 0};


              for (int i = 0; i < Ny -1; i++)
              {
                  for (int j = 0 ; j < Nx -1; j++)
                  {
                      for (int n = 0; n < 2; n++)
                      {
                          int neiy = i + dy[n];
                          int neix = j + dx[n];
                          if(neiy < Ny -1 and neix < Nx -1)
                          {
                              uint8 distance_dahu = compute_dahu( T, depth_node, node_value,  pointxy[i][j],  pointxy[neiy][neix]);
                              std::cout  << i << ", " << j  << " to    "   << neiy  << ","  << neix   << "  dt  " << int(distance_dahu) << std::endl;
                              //std::cout  << pointxy[i][j]  << " to    "   << pointxy[neiy][neix]   << "  dt  " << int(distance_dahu) << std::endl;

                              g.addEdge(i*(Ny-1)+j, neiy*(Ny-1)+neix, distance_dahu);
                          }
                      }
                  }
              }


              int mst_wt = g.kruskalMST();

              g.compute_rnk();

              image2d<uint8> saliency = g.saliencymap();


              io::imsave(saliency, "saliencymap.png");
















//              point2d p1 = point2d(210,160);
//              point2d p2 = point2d(240,160);

////              point2d p3 = point2d(150,137);

////              point2d p4 = point2d(160,60);



//              uint8 distance2 = compute_dahu( T, depth_node, node_value,  p1,  p2);
//              std::cout << int(distance2) << std::endl;

//              image2d<uint8_t> roi_T = roi(T, F, p1*2, p2*2);
//              io::imsave(roi_T, "roi.png");



//              uint8 distance3 = compute_dahu( T, depth_node, node_value,  p1,  p4);
//              std::cout << int(distance3) << std::endl;












//          }




//        }
//    }

//    closedir(directory);

//    std::cout << (float(inter_distance_water_sum)*intra_distance_water.size())/(inter_distance_water.size()*float(intra_distance_water_sum))  << "  "  << (float(inter_distance_dahu_sum)*intra_distance_dahu.size())/(float(inter_distance_dahu.size()*intra_distance_dahu_sum))  << std::endl;

//    std::cout << (float(inter_distance_mst_sum)*intra_distance_mst.size())/(inter_distance_mst.size()*float(intra_distance_mst_sum))  << "  "  << (float(inter_distance_dahu_sum)*intra_distance_dahu.size())/(float(inter_distance_dahu.size()*intra_distance_dahu_sum))  << std::endl;



}

