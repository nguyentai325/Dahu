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







void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb]  \n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}

namespace mln
{

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] != 0  and p[0] != height -1 and p[1] != 0   and p[1]  != width - 1 )
        {
            if (p[0] == 1 or p[0] == height - 2 or p[1] == width -2 or p[1] == 1 )
                return true;
            else
                return false;
        }
        else
                return false;
    }
}

int main(int argc, char** argv)
{
    if (argc < 2)
    usage(argv);
    const char* input_path = argv[1];




    double start_s=clock();


    std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/fastmbd/scr";
    std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/fastmbd/out_fastmbd";
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
            using namespace mln;
            typedef rgb8 V;



            image2d<V> Img;
            io::imread(fileName.c_str(), Img);
            box2d Dom = Img.domain();

            image2d<uint8_t> I(Dom);

            mln_foreach(auto p, Dom)
            {
                I(p) = 0.2989 * Img(p)[0] + 0.5870 * Img(p)[1] + 0.1140 * Img(p)[2];
            }

            //    OUtput : MBD map D
            //    Auxiliaries : U L M Q

            //    State of the pixel
            //    Set M(s) to flooded and initiate for s in S, other pixels are droughty


            // set of seeds S
            unsigned height = I.nrows();
            unsigned width = I.ncols();

            // MBD map, M states
            // U,L


            image2d<uint8_t>  M(Dom);
            image2d<uint8_t> D(Dom);
            image2d<uint8_t> U(Dom);
            image2d<uint8_t> L(Dom);


            for (int i = 0; i < height ; i++)
            {
                for(int j = 0; j < width ; j++)
                {
                    if (i == 0 or i == height - 1 or j == 0 or j == width -1)
                    {
                        M(point2d(i,j)) = 2;
                        D(point2d(i,j)) = 0;
                    }
                    else
                    {
                        M(point2d(i,j)) = 0;
                        D(point2d(i,j)) = 255;
                    }

                    U(point2d(i,j)) = I((point2d(i,j)));
                    L(point2d(i,j)) = I((point2d(i,j)));

                }
            }






            std::vector<std::queue<point2d> > Q(255); // vector of queues

            //    vec[0].push(point2d(1,1)); // push 1 into queue number 0.
            //    std::cout << vec.size()  << std::endl;
            //    point2d point = vec[0].front();
            //    vec[0].pop();
            //    std::cout << point << std::endl;


            //foreach droughty pixel v which has nearest neighbor seed pixel

            int dx[4] = {1 ,-1 , 0 , 0};
            int dy[4] = {0 , 0, 1, -1};


            for (int i = 0; i < height ; i++)
            {
                for(int j = 0; j < width ; j++)
                {
                    if (is_border(point2d(i,j),height,width))
                    {
                        int min_d = 255;
                        int min_n1 = 0;
                        for (int n1 = 0 ; n1 < 4 ; n1++)
                        {
                            int x  = i + dx[n1];
                            int y  = j + dy[n1];
                            if (M(point2d(x,y))==2)
                            {
                                int d = std::abs(I(point2d(i,j)) - I(point2d(x,y)));
                                if (d < min_d)
                                {
                                    min_d = d;
                                    min_n1 = n1;
                                }
                            }
                        }
                        // flow s -> v
                        U(point2d(i,j)) = std::max(U(point2d(i+dx[min_n1],j+dy[min_n1])),I(point2d(i,j)));
                        L(point2d(i,j)) = std::min(L(point2d(i+dx[min_n1],j+dy[min_n1])),I(point2d(i,j)));
                        D(point2d(i,j)) = U(point2d(i,j)) - L(point2d(i,j));
                        // M(v) = waiting
                        M(point2d(i,j)) = 1;

                        // append v to Q[D(v)]
                        Q[D(point2d(i,j))].push(point2d(i,j));
                    }
                }
            }

            //while Q # empty

            int n = 0;
            while (n < 255)
            {
                while (!Q[n].empty())
                {
                    point2d r = Q[n].front();
                    Q[n].pop();

                    if (M(r) == 2)
                        continue;

                    for (int n1 = 0; n1 < 4 ; n1++)
                    {
                        int tx = r[0] + dx[n1];
                        int ty = r[1] + dy[n1];

                        if (M(point2d(tx,ty)) == 0)
                        {
                            U(point2d(tx,ty)) = std::max(U(r),I(point2d(tx,ty)));
                            L(point2d(tx,ty)) = std::min(L(r),I(point2d(tx,ty)));
                            D(point2d(tx,ty)) = U(point2d(tx,ty)) - L(point2d(tx,ty));
                            M(point2d(tx,ty)) = 1;
                            Q[D(point2d(tx,ty))].push(point2d(tx,ty));

                        }



                        else if (M(point2d(tx,ty))==1 and D(point2d(tx,ty))>D(r))
                        {
                            if (D(point2d(tx,ty))> std::max(U(r),I(point2d(tx,ty))) - std::min(L(r),I(point2d(tx,ty))))
                            {
                                U(point2d(tx,ty)) = std::max(U(r),I(point2d(tx,ty)));
                                L(point2d(tx,ty)) = std::min(L(r),I(point2d(tx,ty)));
                                D(point2d(tx,ty)) = U(point2d(tx,ty)) - L(point2d(tx,ty));
                                Q[D(point2d(tx,ty))].push(point2d(tx,ty));

                            }
                        }
                        else
                            continue;
                    }
                }
                n = n +1;
            }

            fileName = outputDirectory + "/" + std::string(_dirent->d_name);

            io::imsave(D, fileName.c_str());

          }
      }

      closedir(directory);

      double stop_s=clock();
      std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;




}
