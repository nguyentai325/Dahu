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

    std::vector<point2d> shortest_path(image2d<uint8_t> grid, point2d src,  point2d dest)
    {
        box2d D = grid.domain();

        unsigned height = grid.nrows();
        unsigned width = grid.ncols();

        //std::cout << height  << "  "  << width << std::endl;



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
        //std::cout << foundDest << std::endl;
        std::vector<point2d>  shortestPath = tracePath (parent, dest);

        return shortestPath;
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


  tbb::task_scheduler_init init;





  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<V> ima;
  io::imread(input_path, ima);



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
    //write_graphviz(file, T);

  }




    unsigned nodenumbers = T.nodes().size();

    auto& nodes1 = T._get_data()->m_nodes;
    auto& S     = T._get_data()->m_S;
    auto& pmap  = T._get_data()->m_pmap;



    box2d D = pmap.domain();


    //std::cout << nodenumbers << std::endl;
    point2d p_topleft = point2d(2,2);
    point2d p_bottomright = point2d(height-100,width-100);
    //std::cout << "bottom right"  << p_bottomright << std::endl;

    point2d p_topright = point2d(2,width-3);
    point2d p_bottomleft = point2d(height -3, 2);



    // compute depth tree

    std::vector<unsigned>  depth_node(nodenumbers);
    std::fill (depth_node.begin(),depth_node.end(),0);


    mln_foreach(auto x, T.nodes())
    {
        if (x.id() != 0)
            depth_node[x.id()] = depth_node[x.get_parent_id()] + 1;
    }


    // Initial step
    // std::cout << "initial  step"   << std::endl;

    std::vector<V>  color_node(nodenumbers);


    mln_foreach(auto x, T.nodes())
    {
        std::vector<V> listpixel1;
        mln_foreach (auto p, x.proper_pset())
            listpixel1.push_back(ima_compo(ima_compo.point_at_index(p)));
        std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end(),lexicographicalorder_less<value_t>());

        V medianvalue1;
        medianvalue1 = listpixel1[listpixel1.size()/2];
        color_node[x.id()] = medianvalue1;
    }




    // find path and lca

    std::vector<unsigned> path;
    V Dahu_tree;
    V min_node_tree;
    V max_node_tree;

    unsigned node_topleft = pmap(p_topleft);
    unsigned node_bottomright = pmap(p_bottomright);
    unsigned node_topright = pmap(p_topright);
    unsigned node_bottomleft = pmap(p_bottomleft);

    unsigned node_start = node_topleft;
    unsigned node_dest = node_bottomright;


    unsigned node_cur = depth_node[node_start] >= depth_node[node_dest] ? node_start : node_dest;
    unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest :node_start  ;



    unsigned d_cur = depth_node[node_cur];
    unsigned d_obj = depth_node[node_obj];


    //std::cout << d_cur   << "  "  << d_obj  << std::endl;

    min_node_tree = color_node[node_cur];
    max_node_tree = color_node[node_cur];

    //std::cout << "ola  "  << std::endl;
    while (d_cur > d_obj)
    {
        d_cur -= 1;
        path.push_back(node_cur);

        //std::cout << int(color_node[node_cur][0]) << "  " << int(color_node[node_cur][1]) << "  " << int(color_node[node_cur][2]) << std::endl;

        for (int k = 0 ; k < 3 ; k++)
        {
            if (color_node[node_cur][k] <  min_node_tree[k])
                min_node_tree[k] = color_node[node_cur][k];
            if (color_node[node_cur][k] >  max_node_tree[k])
                max_node_tree[k] = color_node[node_cur][k];
        }

        node_cur = T.get_node(node_cur).get_parent_id();
    }



    while (node_cur != node_obj)
    {
        path.push_back(node_cur);
        //std::cout << int(color_node[node_cur][0]) << "  " << int(color_node[node_cur][1]) << "  " << int(color_node[node_cur][2]) << std::endl;

        for (int k = 0 ; k < 3 ; k++)
        {
            if (color_node[node_cur][k] <  min_node_tree[k])
                min_node_tree[k] = color_node[node_cur][k];
            if (color_node[node_cur][k] >  max_node_tree[k])
                max_node_tree[k] = color_node[node_cur][k];
        }
        path.push_back(node_obj);
        //std::cout << int(color_node[node_obj][0]) << "  " << int(color_node[node_obj][1]) << "  " << int(color_node[node_obj][2]) << std::endl;

        for (int k = 0 ; k < 3 ; k++)
        {
            if (color_node[node_obj][k] <  min_node_tree[k])
                min_node_tree[k] = color_node[node_obj][k];
            if (color_node[node_obj][k] >  max_node_tree[k])
                max_node_tree[k] = color_node[node_obj][k];
        }
        node_cur = T.get_node(node_cur).get_parent_id();
        node_obj = T.get_node(node_obj).get_parent_id();
    }

    path.push_back(node_cur);
    //std::cout << int(color_node[node_cur][0]) << "  " << int(color_node[node_cur][1]) << "  " << int(color_node[node_cur][2]) << std::endl;

    for (int k = 0 ; k < 3 ; k++)
    {
        if (color_node[node_cur][k] <  min_node_tree[k])
            min_node_tree[k] = color_node[node_cur][k];
        if (color_node[node_cur][k] >  max_node_tree[k])
            max_node_tree[k] = color_node[node_cur][k];
        Dahu_tree[k] = max_node_tree[k] - min_node_tree[k];
    }
    unsigned lca = node_cur;
    
    std::cout << "lca  "  << lca  << std::endl;



    // show tree path

    image2d<uint8_t> path_ima(D);
    extension::fill(path_ima, 0);

    std::vector<uint8_t>  path_node(nodenumbers);
    std::fill (path_node.begin(),path_node.end(),0);

    for(int index = 0 ; index < path.size() ; ++index)
    {
        path_node[path[index]] = 255;
        //std::cout << path[index]  << std::endl;
    }

    mln_foreach(auto x, T.nodes())
    {
		std::cout << x.id()  << std::endl;
        mln_foreach (auto p, x.proper_pset())
        {

            path_ima(F.point_at_index(p)) = path_node[x.id()];

        }
    }


    image2d<uint8_t> under_ima(D);
    extension::fill(under_ima, 0);

    std::vector<uint8_t>  under_node(nodenumbers);
    std::fill (under_node.begin(),under_node.end(),0);


    mln_foreach(auto x, T.nodes())
    {
        if (path_node[x.get_parent_id()] == 255 and path_node[x.id()] == 0)
        {
            under_node[x.id()] = 255;
        }
    }

    mln_foreach(auto x, T.nodes())
    {
        if (under_node[x.get_parent_id()] == 255 and under_node[x.id()] == 0)
        {
            under_node[x.id()] = 255;
        }
        mln_foreach (auto p, x.proper_pset())
        {
            under_ima(F.point_at_index(p)) = under_node[x.id()];
        }
    }
    
    io::imsave(under_ima, "dkm.png");


    image2d<uint8_t> res_ima(D);
    extension::fill(res_ima, 0);

    std::vector<uint8_t>  res_node(nodenumbers);
    std::fill (res_node.begin(),res_node.end(),0);
    res_node[lca] = 255;

    mln_foreach(auto x, T.nodes())
    {
        if (res_node[x.get_parent_id()] == 255 and res_node[x.id()] == 0)
        {
            res_node[x.id()] = 255;
        }
        mln_foreach (auto p, x.proper_pset())
        {
            res_ima(F.point_at_index(p)) = res_node[x.id()];
        }
    }

    //io::imsave(res_ima, "res_ima.png");

    rect2d r = make_rectangle2d(3, 3);
    res_ima = morpho::structural::dilate(res_ima, r);

    io::imsave(res_ima, "vkl.png");


    image2d <uint8_t>  result_ima(D);
    extension::fill(result_ima, 0);

    mln_foreach (auto p , D)
    {
        if (res_ima(p) != under_ima(p))
        {
            result_ima(p) =  255 ;
        }
    }

    io::imsave(result_ima, "result_ima.pgm");


    std::vector<point2d> shortestPath = shortest_path(result_ima, p_topleft, p_bottomright);



    image2d<uint8_t> pathimage(D);


    V truth_dahu_distance;
    V min_path = ima_compo(shortestPath[0]);
    V max_path = ima_compo(shortestPath[0]);

    for (int j = 0 ; j < shortestPath.size() ; j++)
    {
        point2d p = shortestPath[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage(p) = 255;
        for (int k = 0 ; k < 3 ; k++)
        {
            if (ima_compo(p)[k] <  min_path[k])
                min_path[k] = ima_compo(p)[k];
            if (ima_compo(p)[k] > max_path[k])
                max_path[k] = ima_compo(p)[k];
            truth_dahu_distance[k] = max_path[k] - min_path[k];
        }
    }


    //std::cout << int(truth_dahu_distance[0])  << "   "<<int(truth_dahu_distance[1]) <<"  "<<  int(truth_dahu_distance[2]) << std::endl;

    //std::cout << int(Dahu_tree[0])  << "   "<<int(Dahu_tree[1]) <<"  "<<  int(Dahu_tree[2]) << std::endl;

    io::imsave(pathimage, "path.pgm");


}
