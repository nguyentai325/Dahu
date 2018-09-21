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
#include "function.hpp"



//// Number of vertices in the graph
//#define V

namespace mln
{


//    bool is_border (point2d p, unsigned height , unsigned width)
//    {
//        if (p[0] == 0 or p[0]== height -1 or p[1] == 0 or p[1] == width-1)
//            return true;
//        else
//            return false;

//    }

    image2d<rgb8> addborder_color(const image2d<rgb8> ima)
    {
        image2d<rgb8> out(ima.nrows() + 2, ima.ncols() + 2);

        {
          box2d box = ima.domain();
          box.pmin += 1; box.pmax += 1;
          copy(ima, out | box);
        }

        rgb8 median;
        unsigned ncols = ima.ncols(), nrows = ima.nrows();
        {
          std::vector<rgb8> border;
          border.reserve(2 * (nrows + ncols) - 4);
          for (unsigned i = 0; i < ncols; ++i)
          {
            border.push_back(ima.at(0,i));
            border.push_back(ima.at(nrows-1,i));
          }

          for (unsigned i = 1; i < nrows-1; ++i)
          {
            border.push_back(ima.at(i,0));
            border.push_back(ima.at(i,ncols-1));
          }

        std::partial_sort(border.begin(), border.begin() + border.size()/2+1, border.end(), lexicographicalorder_less<rgb8>());//
        if (border.size() % 2 == 0)
            median = border[border.size()/2];
        else
            median = border[border.size()/2];
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

    // A utility function to print the constructed MST stored in parent[]
    void printMST(image2d<point2d> &parent)
    {
       mln_foreach(auto p, parent.domain())
          std::cout <<  p <<   "  " << parent(p) << std::endl;
    }

    // Function to construct and print MST for a graph represented using adjacency
    // matrix representation
//    image2d<point2d> primMST_color(image2d<rgb8> I, std::vector<point2d> &S)
//    {
//         image2d<point2d> parent(I.domain()); // Array to store constructed MST
//         image2d<int> key(I.domain());   // Key values used to pick minimum weight edge in cut
//         image2d<bool> mstSet(I.domain()); // To represent set of vertices not yet included in MST
//         image2d<int> state(I.domain());

//         unsigned height = I.nrows();
//         unsigned width = I.ncols();


//        mln_foreach(auto p, I.domain())
//        {
//             key(p) = INT_MAX;
//             mstSet(p) = false;
//             state(p) = 0;
//         }

//        typedef std::pair<int, point2d> iPair;

//        //std::queue<point2d>  Q;
//        std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > Q;

//        point2d p = point2d(0,0);
//        key(p) = 0;
//        parent(p) = p;
//        Q.push(std::make_pair(0,p));


//        int x[4] = {-1,1,0,0};
//        int y[4] = {0,0,-1,1};



//        while (!Q.empty())
//        {
//            point2d u = Q.top().second;
//            Q.pop();

//            if (mstSet(u) == true)
//                    continue;

//            mstSet(u) = true;
//            S.push_back(u);

//            for (int t = 0; t < 4; t++)
//            {
//                 int x_n = u[0] + x[t];
//                 int y_n = u[1] + y[t];

//                 if( (x_n >= 0 && x_n < height) && (y_n >= 0 && y_n < width) )
//                 {
//                     point2d v = point2d(x_n,y_n);
//                     int temp = std::abs(int(I(v)[0])-int(I(u)[0])) + std::abs(int(I(v)[1])-int(I(u)[1])) + std::abs(int(I(v)[2])-int(I(u)[2]));

//                     if (mstSet(v) == false and temp < key(v))
//                     {
//                         parent(v) = u;
//                         key(v) = temp;
//                         Q.push(std::make_pair(key(v),v));
//                     }

//                 }
//            }
//        }


//        return parent;
//    }
}


// driver program to test above function
int main(int argc, char** argv)
{


    using namespace mln;

    std::string inputDirectory = "/media/minh/DATA/Study/database/MSRA10K_Imgs_GT/MSRA10K_Imgs_GT/Imgs_resize";
    std::string outputDirectory = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/MSTMBD/color";
    std::string outputDirectory_scalar = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/MSTMBD/scalar";

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
            //std::string fileName = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_image/dog_fil.png";

            //const char* input_path = argv[1];

            std::string name = std::string(_dirent->d_name);

            std::istringstream iss(name);
            std::vector<std::string> tokens;
            std::string token;
            while (std::getline(iss, token, '.')) {
                if (!token.empty())
                    tokens.push_back(token);
            }

            std::cout << tokens[1]  << std::endl;
            std::string surname = "jpg";

            if (surname.compare(tokens[1]) == 0)
            {


            image2d<rgb8> ima;
            io::imread(fileName.c_str(), ima);

            image2d<rgb8> a = addborder_color(ima);



            int height = a.nrows();
            int width = a.ncols();

            //// Minimum spanning tree ///////////////////

            std::cout << "minimum spanning tree  "  << std::endl;
            std::vector<point2d> S;
            image2d<point2d> parent = primMST_color(a, S, point2d(0,0));
        //    printMST(parent);

        //    for(int i = 0; i < S.size(); i++)
        //        std::cout << S[i]  << std::endl;


            //// Compute MBD on the MST  /////////////////

            // initiation

            image2d<bool> dejavu(a.domain());
            image2d<bool> seed_node(a.domain());
            image2d<rgb8> min_node(a.domain());
            image2d<rgb8> max_node(a.domain());
            image2d<rgb8> min_temp_node(a.domain());
            image2d<rgb8> max_temp_node(a.domain());
            image2d<unsigned> dist_gray_temp_node(a.domain());
            // distance map
            image2d<rgb8> distance_map(a.domain());
            image2d<uint8> distance_map_scalar(a.domain());


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
                rgb8 min_temp = min_node(parent(S[i]));
                rgb8 max_temp = max_node(parent(S[i]));

                if (dejavu(S[i]) == true)
                {
                    point2d par_p = parent(S[i]);

                    if (seed_node(par_p) == false)
                    {

                        for (int k = 0; k< 3; ++k)
                        {
                            if (min_node(par_p)[k] > min_temp_node(S[i])[k])
                                min_temp[k] = min_temp_node(S[i])[k];

                            if (max_node(par_p)[k] < max_temp_node(S[i])[k])
                                max_temp[k] = max_temp_node(S[i])[k];

                            dis_temp = dis_temp + max_temp[k]  - min_temp[k];
                        }



                        if (dis_temp < dist_gray_temp_node(par_p))
                        {
                            dist_gray_temp_node(par_p) = dis_temp;
                            min_temp_node(par_p) = min_temp;
                            max_temp_node(par_p) = max_temp;
                        }
                        dejavu(par_p) =true;
                        //std::cout << S[i] << "   " << par_p  << "   "<< dis_temp << std::endl;
                    }

                }
            }


            // top down

            std::cout << "chay down xuong "  << std::endl;


            for(int i = 0; i < height*width ; i++)
            {
                rgb8 min_temp = min_node(S[i]);
                rgb8 max_temp = max_node(S[i]);
                unsigned dis_temp = 0;

                if(seed_node(S[i]) == false)
                {

                    for (int k = 0; k< 3; ++k)
                    {
                        if (min_node(S[i])[k] > min_temp_node(parent(S[i]))[k])
                            min_temp[k] = min_temp_node(parent(S[i]))[k];
                        if (max_node(S[i])[k] < max_temp_node(parent(S[i]))[k])
                            max_temp[k] = max_temp_node(parent(S[i]))[k];
                        dis_temp = dis_temp + max_temp[k]  - min_temp[k];
                    }


                    if (dis_temp < dist_gray_temp_node(S[i]))
                    {
                        dist_gray_temp_node(S[i]) = dis_temp;
                        max_temp_node(S[i]) = max_temp;
                        min_temp_node(S[i]) = min_temp;
                    }
                }
            }


            unsigned alpha = 30;
            unsigned beta = 220;



            mln_foreach(auto p, a.domain())
            {
                for (int k = 0; k< 3; ++k)
                {
                    distance_map(p)[k] = max_temp_node(p)[k] - min_temp_node(p)[k];
                }
                distance_map_scalar(p) = distance_map(p)[0]/3 + distance_map(p)[1]/3 + distance_map(p)[2]/3;
            }


//            int max = 0;
//            mln_foreach(auto p, a.domain())
//            {
//                if (distance_map_scalar(p) > max)
//                    max = distance_map_scalar(p);
//            }

//            mln_foreach(auto p, a.domain())
//            {
//                distance_map_scalar(p) = uint8(float(distance_map_scalar(p))/ float(max) *255);
//            }

//            mln_foreach(auto p, a.domain())
//            {
//                if (distance_map_scalar(p) <= alpha)
//                    distance_map_scalar(p) = 0;
//                else if (distance_map_scalar(p) >= beta)
//                    distance_map_scalar(p) = 255;
//                else
//                {
//                    distance_map_scalar(p) = uint8(float(distance_map_scalar(p) - alpha)/ float(beta - alpha) *255);
//                }
//            }


            fileName = outputDirectory + "/" + std::string(_dirent->d_name);
            io::imsave(distance_map, fileName.c_str());

            fileName = outputDirectory_scalar + "/" + std::string(_dirent->d_name);
            io::imsave(distance_map_scalar, fileName.c_str());
//              image2d<rgb8> mst = mbd_mst_distance(a);
//                io::imsave(mst, "map_mst.png");
            }
        }
    }

    closedir(directory);





    return 0;
}
