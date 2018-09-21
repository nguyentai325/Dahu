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
#include <mln/colors/lab.hpp>



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

    template <typename T>
    using lab = internal::vec_base<T, 3, lab_tag>;
    typedef lab<float> lab_f;

    class Lab
    {
        public:

        Lab() {}
        Lab(float L, float a, float b) : L_(L), a_(a), b_(b) {}
        float  L() const { return L_; }
        float& L() { return L_; }
        float  a() const { return a_; }
        float& a() { return a_; }
        float  b() const { return b_; }
        float& b() { return b_; }
        private:
        float L_, a_, b_;
    };

    void rgb8_to_lab(const rgb8& c, float& L, float& a, float& b)
    {
        float
        R = float(c[0])   / 255.0,
        G = float(c[1]) / 255.0,
        B = float(c[2])  / 255.0;

        float X, Y,Z;
        float r1,g1,b1;
        float epsilon = 0.008856;	//actual CIE standard
        float kappa   = 903.3;		//actual CIE standard
        float Xr = 0.950456;	//reference white
        float Yr = 1.0;		//reference white
        float Zr = 1.088754;	//reference white
        double xr,yr,zr;
        double fx, fy, fz;

        if(R <= 0.04045)	r1 = R/12.92;
        else				r1 = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g1 = G/12.92;
        else				g1 = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b1 = B/12.92;
        else				b1 = pow((B+0.055)/1.055,2.4);


        X = r1*0.4124564 + g1*0.3575761 + b1*0.1804375;
        Y = r1*0.2126729 + g1*0.7151522 + b1*0.0721750;
        Z = r1*0.0193339 + g1*0.1191920 + b1*0.9503041;

        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;


        L = 116.0*fy-16.0;
        a = 500.0*(fx-fy);
        b = 200.0*(fy-fz);

    }

    void
    split_Lab(const image2d<Lab>& input,
            image2d<float>& L, image2d<float>& a, image2d<float>& b)
    {
        box2d D = input.domain();
        image2d<float> L_(D), a_(D), b_(D);

        mln_foreach(auto p, D)
        {
            L_(p) = input(p).L();
            a_(p) = input(p).a();
            b_(p) = input(p).b();
        }
        L = L_;
        a = a_;
        b = b_;
    }



    image2d<Lab>
    merge_Lab(const image2d<float>& L, const image2d<float>& a, const image2d<float>& b)
    {
        box2d D = L.domain();
        image2d<Lab> output(D);

        mln_foreach(auto p, D)
        {
            output(p).L() = L(p);
            output(p).a() = a(p);
            output(p).b() = b(p);
        }
        return output;
    }

    image2d<Lab> convert_and_shrink(const image2d<rgb8>& input)
    {
        image2d<Lab> output(input.nrows() ,
                            input.ncols() );
        box2d dom = output.domain();


        mln_foreach(auto p,dom)
        {
            float min = 1e10, max = -1e10, sum = 0;

            float L, a, b;

            rgb8_to_lab(input(p), L, a, b);
            if (L > max)
              max = L;
            if (a < min)
              min = a;
            sum += b;

            output(p).L() = max;
            output(p).a() = min;
            output(p).b() = sum;
        }

        return output;
    }


//    image2d<lab_f> convert_rgb_2_lab( const image2d<rgb8>& input)
//    {
//        image2d<lab_f> output(input.nrows() ,
//                            input.ncols() );
//        box2d dom = output.domain();


//        mln_foreach(auto p,dom)
//        {


//            output(p)  = rgb2lab(input(p));

//        }

//        return output;
//    }


//    image2d<rgb<float> > convert_lab_2_rgb(const image2d<lab_f>& input_lab)
//    {
//        image2d<rgb<float> > output(input_lab.nrows() ,
//                            input_lab.ncols() );
//        box2d dom = output.domain();


//        mln_foreach(auto p,dom)
//        {
//            output(p)  = lab2rgb(input_lab(p));
//            std::cout << output(p)[0] << std::endl;

//        }

//        return output;
//    }


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

    // ************************* convert rgb to Lab  ******************************

    typedef lab<float> lab_f;
//    image2d<lab_f> input_Lab = convert_rgb_2_lab(ima);

//    std::cout << input_Lab(point2d(0,0))  << std::endl;
//    image2d<rgb<float> > asd = convert_lab_2_rgb(input_Lab);
    V color = {255,255,255};
    lab_f out   = rgb2lab(color);
    rgb<float> out1 = lab2rgb(out);

    std::cout << out << std::endl;
    std::cout << out1 << std::endl;

//    io::imsave(asd, "vkl.png");

//    std::cout << (point2d(0,0))  << std::endl;

//    image2d<lab_f> f = addborder(input_Lab, std::less<lab_f>());

//    image2d<lab_f> F = interpolate_k1(f);
//    unsigned height = F.nrows();
//    unsigned width = F.ncols();

//    image2d<lab_f> ima_compo(F.domain() );


//    // COmpute the color tree of shapes T
//    tree_t T;
//    {
//      // 1. Compute the marginal ToS and filter them if necessary.
//      tree_t t[NTREE];
//      tbb::parallel_for(0, (int)NTREE, [&t,&f,a0](int i){
//          t[i] = morpho::cToS(imtransform(f, [i](lab_f x) { return x[i]; }), c4);
//          if (a0 > 0) {
//            grain_filter_inplace(t[i], a0);
//            t[i].shrink_to_fit();
//          }
//        });


//      auto& U0  = t[0]._get_data()->m_Uv;
//      auto& U1  = t[1]._get_data()->m_Uv;
//      auto& U2  = t[2]._get_data()->m_Uv;


//      mln_foreach(point2d p1, U0.domain())
//      {
//          ima_compo(p1)[0] = U0(p1);
//          ima_compo(p1)[1] = U1(p1);
//          ima_compo(p1)[2] = U2(p1);
//      }

//      // 2. Compute the Gos.
//      MyGraph g2;
//      std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
//      std::tie(g2, tlink) = compute_g2<NTREE>(t);

//      // 3. Compute the depth image
//      boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
//      image2d<uint16> imdepth = imchvalue<uint16>(F);

//      write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

//      // debug
//      // io::imsave(imdepth, "depth.tiff");

//      // 4. Compute the saturated maxtree
//      std::tie(T, std::ignore) = satmaxtree(imdepth);
//      T = tree_keep_2F(T);
//      T._reorder_pset();





//      std::ofstream file("abc.dot",
//               std::ios_base::out | std::ios_base::trunc);
//      //write_graphviz(file, T);

//    }


//    // Initiation

//    unsigned nodenumbers = T.nodes().size();

//    auto& nodes1 = T._get_data()->m_nodes;
//    auto& S     = T._get_data()->m_S;
//    auto& pmap  = T._get_data()->m_pmap;

//    box2d D = pmap.domain();

//    image2d <lab_f>  distance_map(D);
//    extension::fill(distance_map, 0);



//    std::vector<lab_f>  min_node(nodenumbers);
//    std::vector<lab_f>  max_node(nodenumbers);

//    std::vector<lab_f>  min_temp_node(nodenumbers);
//    std::vector<lab_f>  max_temp_node(nodenumbers);

//    std::vector<lab_f>  dist(nodenumbers);
//    std::vector<float>  dist_gray_temp_node(nodenumbers);


//    std::vector< std::vector<unsigned> > children(nodenumbers);

//    std::vector<int>  seed_node(nodenumbers);
//    std::fill (seed_node.begin(),seed_node.end(),0);
//    seed_node[0] = 1;
//    std::vector<unsigned>  dejavu(nodenumbers);
//    std::fill (dejavu.begin(),dejavu.end(),0);


//    point2d seedpoints;
//    seedpoints = point2d (69,372);
//    //point2d seedpoints1 = point2d (0, 0);
//    //point2d seedpoints2 = point2d (278, 258);

//    //dejavu[pmap(seedpoints)] = 1;
//    //seed_node[pmap(seedpoints)] = 1;
//    //dejavu[pmap(seedpoints1)] = 1;
//    //seed_node[pmap(seedpoints1)] = 1;

//    //dejavu[pmap(seedpoints2)] = 1;
//    //seed_node[pmap(seedpoints2)] = 1;

//    std::cout << "seed "  << pmap(seedpoints)  << std::endl;


//    mln_foreach(auto p, D)
//    {
//        if (is_border(p,height,width))
//        {
//            dejavu[pmap(p)] = 1;
//            seed_node[pmap(p)] = 1;
//        }
//    }




//    mln_foreach(auto x, T.nodes())
//    {
//        unsigned q = x.get_parent_id();
//        children[q].push_back(x.id());

//        unsigned p = x.first_point();


//        min_node[x.id()] = ima_compo(ima_compo.point_at_index(p));
//        min_temp_node[x.id()] = ima_compo(ima_compo.point_at_index(p));
//        max_node[x.id()] = ima_compo(ima_compo.point_at_index(p));
//        max_temp_node[x.id()] = ima_compo(ima_compo.point_at_index(p));
//        dist_gray_temp_node[x.id()] = 255 *3;


//    }

//    // chay up len


//    std::cout << "chay up len  "  << std::endl;

//    mln_reverse_foreach(auto x, T.nodes())
//    {
//        float dis_temp = 0;
//        lab_f min_temp = min_node[x.get_parent_id()];
//        lab_f max_temp = max_node[x.get_parent_id()];

//        if (dejavu[x.id()] == 1)
//        {
//            unsigned par_p = x.get_parent_id();

//            if (seed_node[par_p] == 0)
//            {

//                for (int i = 0; i< 3; ++i)
//                {
//                    if (min_node[par_p][i] > min_temp_node[x.id()][i])
//                        min_temp[i] = min_temp_node[x.id()][i];

//                    if (max_node[par_p][i] < max_temp_node[x.id()][i])
//                        max_temp[i] = max_temp_node[x.id()][i];

//                    dis_temp = dis_temp + max_temp[i]  - min_temp[i];

//                }

//                if (dis_temp < dist_gray_temp_node[par_p])
//                {
//                    dist_gray_temp_node[par_p] = dis_temp;
//                    min_temp_node[par_p] = min_temp;
//                    max_temp_node[par_p] = max_temp;
//                }
//                dejavu[par_p] =1;
//                //std::cout << x.id() << "   " << par_p  << "   "<< dis_temp << std::endl;

//            }
//        }
//    }

//    // chay down xuong


//    std::cout << "chay down xuong "  << std::endl;

//    std::cout << "root " << T.get_root_id()  << std::endl;

//    mln_foreach (auto x, T.nodes())
//    {
//        lab_f min_temp = min_node[x.id()];
//        lab_f max_temp = max_node[x.id()];
//        unsigned dis_temp = 0;

//        if(seed_node[x.id()] == 0)
//        {
//            for (int i = 0; i < 3; ++i)
//            {
//                if (min_node[x.id()][i] > min_temp_node[x.get_parent_id()][i])
//                    min_temp[i] = min_temp_node[x.get_parent_id()][i];
//                if (max_node[x.id()][i] < max_temp_node[x.get_parent_id()][i])
//                    max_temp[i] = max_temp_node[x.get_parent_id()][i];
//                dis_temp = dis_temp + max_temp[i]  - min_temp[i];
//            }

//            if (dis_temp < dist_gray_temp_node[x.id()])
//            {
//                dist_gray_temp_node[x.id()] = dis_temp;
//                max_temp_node[x.id()] = max_temp;
//                min_temp_node[x.id()] = min_temp;
//                std::cout << max_temp_node[x.id()]- min_temp_node[x.id()]  << std::endl;
//            }
//        }
//    }


//    // back propagation

//    mln_foreach(auto x, T.nodes())
//    {
//        mln_foreach (auto p, x.proper_pset())
//        {
//            distance_map(F.point_at_index(p)) = max_temp_node[x.id()] - min_temp_node[x.id()];
//        }
//    }

//    image2d<rgb8> dahu_distance = convert_lab_2_rgb(distance_map);

//    io::imsave(dahu_distance, "dmap1.png");









}
