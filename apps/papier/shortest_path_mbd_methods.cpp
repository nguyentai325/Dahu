#include<bits/stdc++.h>


#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
# include <mln/core/extension/fill.hpp>


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

# include <mln/morpho/structural/dilate.hpp>


#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"

#include <ctime>
#include <math.h>
#include <float.h>
#include <time.h>
#include <chrono>


#include <mln/morpho/tos/ctos.hpp>


#include <functional>
#include <queue>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "function.hpp"



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
    bool compare (int i,int j) { return (i<j); }


    bool is_top (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 )
            return true;
        else
            return false;

    }

//    bool is_border (point2d p, unsigned height , unsigned width)
//    {
//        if (p[0] != 0  and p[0]  != 1  and p[0] != height -1 and p[0] != height -2 and p[1] != 0  and p[1]  != 1 and p[1]  != width - 1 and p[1] != width -2)
//        {
//            if (p[0] == 2 or p[0] == height - 3 or p[1] == width -3 or p[1] == 2  )
//                return true;
//            else
//                return false;
//        }

//        //if (p[0] == 0 or p[0] == height - 1 or p[1] == width -1 or p[1] == 0  )
//        //	return true;
//        //else
//        //	return false;

//    }

    typedef std::pair<float, point2d > pPair;

    // A Utility Function to check whether given cell (row, col)
    // is a valid cell or not.
    bool isValid(point2d p, box2d D )
    {
        // Returns true if row number and column number
        // is in range
        return (p[0] >= 0) && (p[0] < D.pmax[0]) &&
                (p[1] >= 0) && (p[1] < D.pmax[1]);
    }

    // A Utility Function to check whether the given cell is
    // blocked or not
    bool isUnBlocked(image2d<uint8_t> grid, point2d p)
    {
        // Returns true if the cell is not blocked else false
        if (grid(p) == 255)
            return (true);
        else
            return (false);
    }

    // A Utility Function to check whether destination cell has
    // been reached or not
    bool isDestination(point2d p, point2d dest)
    {
        if ( p == dest)
            return (true);
        else
            return (false);
    }


    // A Utility Function to calculate the 'h' heuristics.
    double calculateHValue(point2d p, point2d dest)
    {
        // Return using the distance formula
        return ((double)sqrt ((p[0]-dest[0])*(p[0]-dest[0])
                              + (p[1]-dest[1])*(p[1]-dest[1])));
    }


    std::vector<point2d> tracePath(image2d<point2d> parent, point2d dest)
    {
        //printf ("\nThe Path is ");
//        int row = dest.first;
//        int col = dest.second;
        point2d p = dest;

        std::vector<point2d> Path;

        while (!(parent(p) == p ))
        {
            Path.push_back(p);
            point2d temp_p = parent(p);
            p = temp_p;
        }

        Path.push_back(p);
//        while (!Path.empty())
//        {
//            point2d p = Path.top();
//            Path.pop();
//            printf("-> (%d,%d) ",p[0],p[1]);
//        }
        return Path;

    }

    std::vector<point2d> tracePath_mst(image2d<point2d> parent, point2d p_bottomright)
    {
        printf ("\nThe Path is ");
//        int row = dest.first;
//        int col = dest.second;
        point2d p = p_bottomright;

        std::vector<point2d> Path;

        while (!(parent(p) == p ))
        {
            Path.push_back(p);
            point2d temp_p = parent(p);
            p = temp_p;
        }

        Path.push_back(p);
//        while (!Path.empty())
//        {
//            point2d p = Path.top();
//            Path.pop();
//            printf("-> (%d,%d) ",p[0],p[1]);
//        }
        return Path;

    }

    std::vector<point2d> shortest_path(image2d<uint8_t> grid, point2d src,  point2d dest)
    {
        box2d D = grid.domain();


        unsigned height = grid.nrows();
        unsigned width = grid.ncols();

        std::cout << height  << "  "  << width << std::endl;



        if (isValid (src, D) == false)
        {
            std::printf ("Source is invalid\n");
        }

        // If the destination is out of range
        if (isValid (dest, D) == false)
        {
            std::printf ("Destination is invalid\n");
        }

        // Either the source or the destination is blocked
        if (isUnBlocked(grid, src) == false ||
                isUnBlocked(grid, dest) == false)
        {
            std::printf ("Source or the destination is blocked\n");
        }

        // If the destination cell is the same as source cell
        if (isDestination(src, dest) == true)
        {
            std::printf ("We are already at the destination\n");
        }


        // Create a closed list and initialise it to false which means
        // that no cell has been included yet
        // This closed list is implemented as a boolean 2D array

        image2d<bool>  closedList(D);

        // Declare a 2D array of structure to hold the details
        //of that cell  : f,g,h , parent

        image2d<float> f(D);
        image2d<float> g(D);
        image2d<float> h(D);
        image2d<point2d> parent(D);

        mln_foreach(auto p, D)
        {
            closedList(p) = false;
            f(p) = FLT_MAX;
            g(p) = FLT_MAX;
            h(p) = FLT_MAX;
            parent(p) = point2d(-1,-1);
        }


        // Initialising the parameters of the starting node

        f(src) = 0.0;
        g(src) = 0.0;
        h(src) = 0.0;
        parent(src) = src;


        /*
             Create an open list having information as-
             <f, <i, j>>
             where f = g + h,
             and i, j are the row and column index of that cell
             Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
             This open list is implenented as a set of pair of pair.*/

        std::set<pPair> openList;
        openList.insert(std::make_pair(0.0,src));

        // We set this boolean value as false as initially
        // the destination is not reached.

        bool foundDest = false;

        while(!openList.empty() and foundDest == false)
        {
            pPair pair = *openList.begin();
            point2d p = pair.second;
            //std::cout << p << std::endl;
            // Remove this vertex from the open list
            openList.erase(openList.begin());

            // Add this vertex to the open list
            closedList(p) = true;

            /*
                Generating all the 8 successor of this cell

                    N.W   N   N.E
                      \   |   /
                       \  |  /
                    W----Cell----E
                         / | \
                       /   |  \
                    S.W    S   S.E

                Cell-->Popped Cell (i, j)
                N -->  North       (i-1, j)
                S -->  South       (i+1, j)
                E -->  East        (i, j+1)
                W -->  West           (i, j-1)
                N.E--> North-East  (i-1, j+1)
                N.W--> North-West  (i-1, j-1)
                S.E--> South-East  (i+1, j+1)
                S.W--> South-West  (i+1, j-1)*/

            // To store the 'g', 'h' and 'f' of the 8 successors
            float gNew, hNew, fNew;

            //----------- 1st Successor (North) ------------
            int dx[8] = {-1, +1,  0 ,  0, -1, -1, +1, +1};
            int dy[8] = { 0,  0, +1 , -1, +1, -1, +1, -1};


            for (int n1 = 0 ; n1 < 8 ; n1++)
            {
                point2d p_new = point2d(p[0]+dx[n1], p[1]+dy[n1]);
                if (isValid(p_new, D) == true)
                {
                    // If the destination cell is the same as the
                    // current successor
                    if (isDestination(p_new, dest) == true)
                    {
                        // Set the Parent of the destination cell
                        parent(p_new) = p;
                        //std::printf ("The destination cell is found\n");
                        //tracePath (parent, dest);
                        foundDest = true;
                        break;
                    }
                    // If the successor is already on the closed
                    // list or if it is blocked, then ignore it.
                    // Else do the following
                    else if (closedList(p_new) == false &&
                             isUnBlocked(grid,p_new) == true)
                    {
                        gNew = g(p) + 1.0;
                        hNew = calculateHValue (p_new, dest);
                        fNew = gNew + hNew;

                        // If it isn’t on the open list, add it to
                        // the open list. Make the current square
                        // the parent of this square. Record the
                        // f, g, and h costs of the square cell
                        //                OR
                        // If it is on the open list already, check
                        // to see if this path to that square is better,
                        // using 'f' cost as the measure.
                        if (f(p_new) == FLT_MAX || f(p_new) > fNew)
                        {
                            openList.insert(std::make_pair(fNew,p_new) );

                            // Update the details of this cell
                            f(p_new) = fNew;
                            g(p_new) = gNew;
                            h(p_new) = hNew;
                            parent(p_new) = p;
                        }
                    }
                }
            }
        }


        // When the destination cell is not found and the open
        // list is empty, then we conclude that we failed to
        // reach the destiantion cell. This may happen when the
        // there is no way to destination cell (due to blockages)
//        if (foundDest == false)
//            std::printf("Failed to find the Destination Cell\n");
//        else
        std::cout << foundDest << std::endl;
        std::vector<point2d>  shortestPath = tracePath (parent, dest);

        return shortestPath;
    }

//    image2d<uint8_t> roi(tree_t T, image2d<rgb8> F, point2d p1, point2d p2)
//    {

//        unsigned nodenumbers = T.nodes().size();

//        auto& nodes1 = T._get_data()->m_nodes;
//        auto& S     = T._get_data()->m_S;
//        auto& pmap  = T._get_data()->m_pmap;

//        box2d D = pmap.domain();

//        std::vector<unsigned>  depth_node(nodenumbers);
//        std::fill (depth_node.begin(),depth_node.end(),0);


//        mln_foreach(auto x, T.nodes())
//        {
//            if (x.id() != 0)
//                depth_node[x.id()] = depth_node[x.get_parent_id()] + 1;
//        }


//        // find path and lca

//        std::vector<unsigned> path;

//        unsigned node_1 = pmap(p1);
//        unsigned node_2 = pmap(p2);

//        unsigned node_start = node_1;
//        unsigned node_dest = node_2;

//        unsigned node_cur = depth_node[node_start] >= depth_node[node_dest] ? node_start : node_dest;
//        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest :node_start  ;

//        unsigned d_cur = depth_node[node_cur];
//        unsigned d_obj = depth_node[node_obj];


//        while (d_cur > d_obj)
//        {
//            d_cur -= 1;
//            path.push_back(node_cur);
//            node_cur = T.get_node(node_cur).get_parent_id();
//        }

//        while (node_cur != node_obj)
//        {
//            path.push_back(node_cur);
//            path.push_back(node_obj);
//            node_cur = T.get_node(node_cur).get_parent_id();
//            node_obj = T.get_node(node_obj).get_parent_id();
//        }

//        path.push_back(node_cur);

//        unsigned lca = node_cur;


//        // show tree path


//        std::vector<uint8_t>  path_node(nodenumbers);
//        std::fill (path_node.begin(),path_node.end(),0);

//        for(int index = 0 ; index < path.size() ; ++index)
//        {
//            path_node[path[index]] = 255;
//        }




//        image2d<uint8_t> under_ima(D);
//        extension::fill(under_ima, 0);

//        std::vector<uint8_t>  under_node(nodenumbers);
//        std::fill (under_node.begin(),under_node.end(),0);


//        mln_foreach(auto x, T.nodes())
//        {
//            if (path_node[x.get_parent_id()] == 255 and path_node[x.id()] == 0)
//            {
//                under_node[x.id()] = 255;
//            }
//        }

//        mln_foreach(auto x, T.nodes())
//        {
//            if (under_node[x.get_parent_id()] == 255 and under_node[x.id()] == 0)
//            {
//                under_node[x.id()] = 255;
//            }
//            mln_foreach (auto p, x.proper_pset())
//            {
//                under_ima(F.point_at_index(p)) = under_node[x.id()];
//            }
//        }


//        image2d<uint8_t> res_ima(D);
//        extension::fill(res_ima, 0);

//        std::vector<uint8_t>  res_node(nodenumbers);
//        std::fill (res_node.begin(),res_node.end(),0);
//        res_node[lca] = 255;

//        mln_foreach(auto x, T.nodes())
//        {
//            if (res_node[x.get_parent_id()] == 255 and res_node[x.id()] == 0)
//            {
//                res_node[x.id()] = 255;
//            }
//            mln_foreach (auto p, x.proper_pset())
//            {
//                res_ima(F.point_at_index(p)) = res_node[x.id()];
//            }
//        }


//        rect2d r = make_rectangle2d(3, 3);
//        res_ima = morpho::structural::dilate(res_ima, r);


//        image2d <uint8_t>  result_ima(D);
//        extension::fill(result_ima, 0);

//        mln_foreach (auto p , D)
//        {
//            if (res_ima(p) != under_ima(p))
//            {
//                result_ima(p) =  255 ;
//            }
//        }


//        return result_ima;
//    }



    template <typename V>
    image2d<V> waterflow_mbd_distance_2_points(image2d<V> I, point2d p_topleft, point2d p_bottomright, V & waterflow_distance, image2d<point2d> & parent )
    {
        box2d Dom = I.domain();

        // set of seeds S
        unsigned height = I.nrows();
        unsigned width = I.ncols();

        image2d<unsigned>  state(Dom);
        typedef std::pair<rgb8,rgb8> pair_t;
        image2d<pair_t> mm(Dom);
        image2d<rgb8> dmap(Dom);
        image2d<uint8> dmap_scalar(Dom);

        //image2d<uint8_t> dmap1(Dom);

        // priority queue

        std::vector<std::queue<point2d> > Q(256*3);




        // put seed on the border of the image
        // change the state of the pixel
        for (int i = 0; i < height ; i++)
        {
            for(int j = 0; j < width ; j++)
            {
                point2d p = point2d(i,j);
                if (p == p_topleft)
                {
                    state(p) = 1;
                    dmap(p) = {0,0,0};
                    Q[dmap(p)[0]+dmap(p)[1]+dmap(p)[2]].push(p);
                    mm(p) = pair_t(I(p),I(p));
                    parent(p) = p;
                }
                else
                {
                    state(p) = 0;
                    dmap(p) = {255,255,255};
                    mm(p) = pair_t(I(p),I(p));
                    parent(p) = point2d(-1,-1);
                }
            }
        }

    //    stop_s=clock();
    //    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;


        int dx[4] = {1 ,-1 , 0 , 0};
        int dy[4] = {0 , 0, 1, -1};

        // proceed the propagation of the pixel from border to the center of image until all of pixels is passed

        for (int lvl = 0; lvl < 256*3 ; lvl++)
        {
            while (!Q[lvl].empty())
            {
                point2d p = Q[lvl].front();
                Q[lvl].pop();
                //std::cout << "p  " << p   << "  state(p)   "<< state(p) << std::endl;
                if (state(p) == 2)
                    continue;

                state(p) = 2;
                //dmap1(p) = lvl;

                for (int n1 = 0 ; n1 < 4 ; n1++)
                {
                    int x  = p[0] + dx[n1];
                    int y  = p[1] + dy[n1];

                    if (x > 0 && x < height && y > 0 and y < width)
                    {
                        point2d r = point2d(x,y);

                        if (state(r)==1 && dmap(r)[0] + dmap(r)[1] + dmap(r)[2] > dmap(p)[0] + dmap(p)[1] + dmap(p)[2])
                        {

                            mm(r) = mm(p);
                            for (int i = 0; i < 3; i++)
                            {
                                if (I(r)[i] < mm(r).first[i])
                                  mm(r).first[i] = I(r)[i];
                                if (I(r)[i] > mm(r).second[i])
                                  mm(r).second[i] = I(r)[i];
                            }
                            if (dmap(r)[0] + dmap(r)[1] + dmap(r)[2] > mm(r).second[0] - mm(r).first[0] + mm(r).second[1] - mm(r).first[1] + mm(r).second[2] - mm(r).first[2])
                            {

                                dmap(r) = mm(r).second - mm(r).first;
                                Q[dmap(r)[0]+dmap(r)[1]+dmap(r)[2]].push(r);
                                parent(r) = p;
                            }
                        }

                        else if (state(r)==0)
                        {

                            mm(r) = mm(p);
                            for (int i = 0 ; i < 3; i++)
                            {
                                if (I(r)[i] < mm(r).first[i])
                                  mm(r).first[i] = I(r)[i];
                                if (I(r)[i] > mm(r).second[i])
                                  mm(r).second[i] = I(r)[i];
                            }

                            dmap(r) = mm(r).second - mm(r).first;
                            Q[dmap(r)[0]+dmap(r)[1]+dmap(r)[2]].push(r);
                            state(r) =1;
                            parent(r) =p;

                        }
                        else
                            continue;

                    }

                }
            }
        }
        waterflow_distance  = dmap(p_bottomright);

        return dmap;


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
    int a0 = std::atoi(argv[2]);
    int a1 = std::atoi(argv[3]);
    int lambda = std::atoi(argv[4]);
    const char* output_path = argv[5];
    //  const char* output_path_t0 = argv[6];
    //  const char* output_path_t1 = argv[7];
    //  const char* output_path_t2 = argv[8];


    tbb::task_scheduler_init init;

    std::string outputDirectory_Dahu_roi = "/media/minh/DATA/Study/Results/Compare_shortest_path/Dahu/roi";
    std::string outputDirectory_Dahu_path = "/media/minh/DATA/Study/Results/Compare_shortest_path/Dahu/path";
    std::string outputDirectory_waterflow_path = "/media/minh/DATA/Study/Results/Compare_shortest_path/Waterflow";
    std::string outputDirectory_mstmbd_path = "/media/minh/DATA/Study/Results/Compare_shortest_path/MSTMBD";


    std::string fileName = input_path;


    // 1. Compute the individual ToS
    using namespace mln;
    typedef rgb8 V;

    image2d<V> ima;
    io::imread(fileName.c_str(), ima);



    image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
    image2d<V> short_test_ima = interpolate_k1(f);
    box2d d = short_test_ima.domain();
    unsigned height = short_test_ima.nrows();
    unsigned width = short_test_ima.ncols();

    image2d<rgb8> pathimage(f.domain());
    image2d<rgb8> pathimage_mst(f.domain());
    mln_foreach( auto p, f.domain())
    {
        pathimage(p) = f(p);
        pathimage_mst(p) = f(p);
    }


    // COmpute the color tree of shapes T
    tree_t t[NTREE];
    tree_t T;
    {

        // 1. Compute the marginal ToS and filter them if necessary.

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
        image2d<uint16> imdepth = imchvalue<uint16>(short_test_ima);

        write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

        // debug
        // io::imsave(imdepth, "depth.tiff");

        // 4. Compute the saturated maxtree
        std::tie(T, std::ignore) = satmaxtree(imdepth);
        T = tree_keep_2F(T);
        T._reorder_pset();


    }


    t[0]._reorder_pset();
    t[1]._reorder_pset();
    t[2]._reorder_pset();


    point2d p_topleft = point2d(195,209);
    point2d p_bottomright = point2d(164,49);




    for (int i = 0; i < ima.nrows(); i++)
    {
        for(int j = 0; j < ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_bottomright[0] - 4 and i < p_bottomright[0] + 4 and j > p_bottomright[1] - 4 and j < p_bottomright[1]+4) )
            {
                ima(point2d(i,j)) = {255,0,0};
            }

        }
    }


    // Shortest path in Dahu distance case ////////////////////////////////////////////////////////////////


    image2d<uint8_t> roi_T = roi(T, short_test_ima, p_topleft*2, p_bottomright*2);


    std::vector<point2d> path_optimal =  shortest_path( roi_T,  p_topleft*2,   p_bottomright*2);

    //image2d<rgb8> short_test_ima(d);

    for (int i = 0; i < short_test_ima.nrows(); i++)
    {
        for(int j = 0; j < short_test_ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0]*2 - 4 and i < p_topleft[0]*2 + 4 and j > p_topleft[1]*2 - 4 and j < p_topleft[1]*2 + 4)  or (i > p_bottomright[0]*2-4 and i < p_bottomright[0]*2+4 and j > p_bottomright[1]*2-4 and j < p_bottomright[1]*2+4) )
            {
                short_test_ima(point2d(i,j)) = {255,0,0};
            }
        }
    }

    for (unsigned i = 0; i < path_optimal.size() ; i++)
    {
        short_test_ima(path_optimal[i]) = {0,0,255};
    }

    std::string fileName_out = outputDirectory_Dahu_roi + "/" + fileName;
    io::imsave(roi_T, fileName_out.c_str());

    std::string fileName_out_path = outputDirectory_Dahu_path + "/" + fileName;
    io::imsave(short_test_ima, fileName_out_path.c_str());



    // Shortest path in Fast MBD case  ////////////////////////////////////////////////////////////////////////////////////






    // Shortest path in MSTMBD case ////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<point2d> S1;
    image2d<point2d> parent_mst = primMST_color(pathimage_mst, S1, p_topleft);


    image2d<unsigned> depth_node_mst = depth_node_mst_compute(parent_mst, S1 );


    std::vector<point2d>  shortestPath_mst = tracePath_mst (parent_mst, p_bottomright);



    for (int i = 0; i < ima.nrows(); i++)
    {
        for(int j = 0; j < ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_bottomright[0]-4 and i < p_bottomright[0]+4 and j > p_bottomright[1]-4 and j < p_bottomright[1]+4) )
            {
                pathimage_mst(point2d(i,j)) = {255,0,0};
            }

        }
    }

    for (int j = 0 ; j < shortestPath_mst.size() ; j++)
    {
        point2d p = shortestPath_mst[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage_mst(p) = {0,0,255};
    }

    fileName_out_path = outputDirectory_mstmbd_path + "/" + fileName;
    io::imsave(pathimage_mst, fileName_out_path.c_str());


    // Shortest path in waterflow case ////////////////////////////////////////////////////////////////////////////////////////////

    V waterflow_distance;
    box2d Dom = pathimage.domain();
    image2d<point2d > parent(Dom);
    image2d<V> waterflow_map = waterflow_mbd_distance_2_points( pathimage,  p_topleft,  p_bottomright, waterflow_distance, parent );

    std::vector<point2d>  shortestPath = tracePath (parent, p_bottomright);


    for (int i = 0; i < ima.nrows(); i++)
    {
        for(int j = 0; j < ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_bottomright[0]-4 and i < p_bottomright[0]+4 and j > p_bottomright[1]-4 and j < p_bottomright[1]+4) )
            {
                pathimage(point2d(i,j)) = {255,0,0};
            }

        }
    }


    for (int j = 0 ; j < shortestPath.size() ; j++)
    {
        point2d p = shortestPath[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage(p)  = {0,0,255};
    }


    fileName_out_path = outputDirectory_waterflow_path + "/" + fileName;
    io::imsave(pathimage, fileName_out_path.c_str());


}
