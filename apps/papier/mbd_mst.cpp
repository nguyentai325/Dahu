// A C / C++ program for Prim's Minimum Spanning Tree (MST) algorithm.
// The program is for adjacency matrix representation of the graph

#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>



#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/dontcare.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/property_map/function_property_map.hpp>

#include <queue>
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
//#include "function.hpp"



//// Number of vertices in the graph
//#define V

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

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 or p[0]== height -1 or p[1] == 0 or p[1] == width-1)
            return true;
        else
            return false;

    }

//    bool is_border (point2d p, unsigned height , unsigned width)
//    {
//        if (p[0] != 0   and p[0] != height -1  and p[1] != 0  and p[1]  != width - 1 )
//        {
//            if (p[0] == 2 or p[0] == height - 2 or p[1] == width - 2 or p[1] == 1  )
//                return true;
//            else
//                return false;
//        }

//    }
    // A utility function to find the vertex with minimum key value, from
    // the set of vertices not yet included in MST

    point2d minKey(image2d<int> &key, image2d<bool> &mstSet)
    {
        int min = INT_MAX;
        point2d min_point;

        mln_foreach(auto p, key.domain())
        {
            if (mstSet(p) == false and key(p) < min)
            {
                min = key(p);
                min_point = p;
            }
        }
        return min_point;
    }

    // A utility function to print the constructed MST stored in parent[]
    void printMST(image2d<point2d> &parent)
    {
       mln_foreach(auto p, parent.domain())
          std::cout <<  p <<   "  " << parent(p) << std::endl;
    }

    // Function to construct and print MST for a graph represented using adjacency
    // matrix representation
    image2d<point2d> primMST(image2d<uint8> I, std::vector<point2d> &S)
    {
         image2d<point2d> parent(I.domain()); // Array to store constructed MST
         image2d<int> key(I.domain());   // Key values used to pick minimum weight edge in cut
         image2d<bool> mstSet(I.domain()); // To represent set of vertices not yet included in MST
         typedef std::pair<unsigned,unsigned> pair_t;
         image2d<pair_t> mm(I.domain());
         image2d<uint8> dmap(I.domain());

         unsigned height = I.nrows();
         unsigned width = I.ncols();


        mln_foreach(auto p, I.domain())
        {
             key(p) = INT_MAX;
             mstSet(p) = false;
             parent(p) = point2d(-1,-1);
         }

        typedef std::pair<int, point2d> iPair;

        //std::queue<point2d>  Q;
        std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > Q;

        point2d p = point2d(0,0);  // source pixel
        key(p) = 0;
        parent(p) = p;
        Q.push(std::make_pair(0,p));
        mm(p) = pair_t(I(p), I(p));


        int x[4] = {-1,1,0,0};
        int y[4] = {0,0,-1,1};



        while (!Q.empty())
        {
            point2d u = Q.top().second;
            Q.pop();

            if (mstSet(u) == true)
                continue;

            mstSet(u) = true;
            S.push_back(u);

            for (int t = 0; t < 4; t++)
            {
                 int x_n = u[0] + x[t];
                 int y_n = u[1] + y[t];

                 if( (x_n >= 0 && x_n < height) && (y_n >= 0 && y_n < width) )
                 {
                     point2d v = point2d(x_n,y_n);

                     if (mstSet(v) == false and std::abs(int(I(v))-int(I(u))) < key(v))
                     {
                         parent(v) = u;
                         key(v) = std::abs(int(I(v))-int(I(u)));
                         Q.push(std::make_pair(key(v),v));

                         mm(v) = mm(u);
                         if (I(v) < mm(v).first)
                            mm(v).first = I(v);
                         if (I(v) > mm(v).second)
                            mm(v).second = I(v);
                         dmap(v) = mm(v).second - mm(v).first;
                     }

                 }
            }
        }


        io::imsave(dmap, "map_mst_mbd1.pgm");


        return parent;
    }




}


// driver program to test above function
int main(int argc, char** argv)
{
    using namespace mln;

    std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_filter";
    std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/output_mstmbd_gray/test_ima";
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
    //std::string fileName = "/home/minh/Documents/code/code_edwin/build/apps/papier/out-027.jpg";



            //// Load image ///////////////////////


            const char* input_path = argv[1];

            image2d<rgb8> ima;
            io::imread(fileName.c_str(), ima);

            image2d<uint8>  ima_gray = rgb2gray(ima);

            image2d<uint8> a = addborder_gray(ima_gray);


            int height = a.nrows();
            int width = a.ncols();

            //// Minimum spanning tree ///////////////////

            std::cout << "minimum spanning tree  "  << std::endl;
            std::vector<point2d> S;
            image2d<point2d> parent = primMST(a, S);
            std::cout << S.size() << std::endl;
           // printMST(parent);

        //    for(int i = 0; i < S.size(); i++)
        //        std::cout << S[i]  << std::endl;


            //// Compute MBD on the MST  /////////////////

            // initiation

            image2d<bool> dejavu(a.domain());
            image2d<bool> seed_node(a.domain());
            image2d<uint8> min_node(a.domain());
            image2d<uint8> max_node(a.domain());
            image2d<uint8> min_temp_node(a.domain());
            image2d<uint8> max_temp_node(a.domain());
            image2d<uint8> dist_gray_temp_node(a.domain());
            // distance map
            image2d<uint8> distance_map(a.domain());



            mln_foreach(auto p, a.domain())
            {
                if (is_border(p,height,width))
                {
                    dejavu(p) = true;
                    seed_node(p) = true;

                }
                else
                {
                    dejavu(p) = false;
                    seed_node(p) = false;
                }

            }


            mln_foreach(auto p, a.domain())
            {

                min_node(p) = a(p);
                min_temp_node(p) = a(p);
                max_node(p) = a(p);
                max_temp_node(p) = a(p);
                dist_gray_temp_node(p) = INT_MAX;


            }



            // bottom up


            std::cout << "chay up len  "  << std::endl;


            for(int i = height*width -1; i >=0 ; i--)
            {
                unsigned dis_temp = 0;
                uint8 min_temp = min_node(parent(S[i]));
                uint8 max_temp = max_node(parent(S[i]));

                if (dejavu(S[i]) == true)
                {
                    point2d par_p = parent(S[i]);

                    if (seed_node(par_p) == false)
                    {


                        if (min_node(par_p) > min_temp_node(S[i]))
                            min_temp = min_temp_node(S[i]);

                        if (max_node(par_p) < max_temp_node(S[i]))
                            max_temp = max_temp_node(S[i]);

                        //std::cout << int(min_node(par_p))  << "   " << int(min_temp_node(S[i])) <<  "  "  << int(min_temp)<< std::endl;
                        //std::cout << int(max_node(par_p))  << "   " << int(max_temp_node(S[i])) <<  "  "  << int(max_temp)<< std::endl;


                        dis_temp = dis_temp + max_temp  - min_temp;

                        //std::cout << S[i]  << "   " << par_p << "  "<< dis_temp << std::endl;

                        if (dis_temp < dist_gray_temp_node(par_p))
                        {
                            dist_gray_temp_node(par_p) = dis_temp;
                            min_temp_node(par_p) = min_temp;
                            max_temp_node(par_p) = max_temp;
                        }
                    dejavu(par_p) =true;
                    }

                }
            }

        //    mln_foreach(auto p, a.domain())
        //    {
        //        distance_map(p) = max_temp_node(p) - min_temp_node(p);
        //        std::cout <<  p <<  "   "  <<  int (distance_map(p))  << std::endl;
        //    }

            // top down

            std::cout << "chay down xuong "  << std::endl;


            for(int i = 0; i < height*width ; i++)
            {
                uint8 min_temp = min_node(S[i]);
                uint8 max_temp = max_node(S[i]);
                unsigned dis_temp = 0;

                if(seed_node(S[i]) == false)
                {

                    if (min_node(S[i]) > min_temp_node(parent(S[i])))
                        min_temp = min_temp_node(parent(S[i]));
                    if (max_node(S[i]) < max_temp_node(parent(S[i])))
                        max_temp = max_temp_node(parent(S[i]));
                    dis_temp = dis_temp + max_temp  - min_temp;


                    if (dis_temp < dist_gray_temp_node(S[i]))
                    {
                        dist_gray_temp_node(S[i]) = dis_temp;
                        max_temp_node(S[i]) = max_temp;
                        min_temp_node(S[i]) = min_temp;
                    }
                }
            }



            mln_foreach(auto p, a.domain())
            {
                distance_map(p) = max_temp_node(p) - min_temp_node(p);
        //        std::cout <<  p <<  "   "  <<  int (distance_map(p))  << std::endl;
            }


            fileName = outputDirectory + "/" + std::string(_dirent->d_name);

            io::imsave(distance_map, fileName);




        }
    }

    closedir(directory);



    return 0;
}
