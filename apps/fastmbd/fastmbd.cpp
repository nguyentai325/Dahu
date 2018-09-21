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










void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb]  \n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}



int main(int argc, char** argv)
{
    if (argc < 2)
    usage(argv);
    const char* input_path = argv[1];
    std::string fileName = input_path;

    using namespace mln;
    typedef rgb8 V;

    double start_s=clock();


    image2d<V> Img;
    io::imread(fileName.c_str(), Img);
    box2d Dom = Img.domain();

    image2d<uint8_t> I(Dom);

    mln_foreach(auto p, Dom)
    {
        I(p) = 0.2989 * Img(p)[0] + 0.5870 * Img(p)[1] + 0.1140 * Img(p)[2];
    }

    // number of passes K
    int K = 3;


    // set of seeds S
    unsigned height = I.nrows();
    unsigned width = I.ncols();

    std::cout << height<< std::endl;
    std::cout << width << std::endl;



    // MBD map



    image2d<uint8_t> D(Dom);
    int i , j;

    for ( i= 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            if (i == 0 )
                D(point2d(i,j)) = 0;
            else if (i == height -1)
                D(point2d(i,j)) = 0;

            else if (j == 0 )
                D(point2d(i,j)) = 0;
            else if (j == width -1 )
                D(point2d(i,j)) = 0;
            else
                D(point2d(i,j)) = 255;
        }
    }



    // U,L

    image2d<uint8_t> U(Dom);
    image2d<uint8_t> L(Dom);


    mln_foreach(auto p , Dom)
    {
        U(p) = I(p);
        L(p) = I(p);
    }


    for (int k = 0; k < K; k++)
    {
        if (k % 2 == 0)
        {
            uint8_t ix, u1, l1, b1, u2, l2, b2, d;

            for (i = 1; i < height - 1 ; i++)
            {
                for (j = 1; j < width - 1; j++)
                {
                    point2d p = point2d(i,j);
                    ix = I(p);
                    u1 = U(point2d(i-1,j));
                    l1 = L(point2d(i-1,j));
                    b1 = std::max(u1,ix) - std::min(l1,ix);

                    u2 = U(point2d(i,j-1));
                    l2 = L(point2d(i,j-1));
                    b2 = std::max(u2,ix) - std::min(l2,ix);


                    d = D(p);



                    if (d <= b1 && d <= b2)
                        continue;
                    else if (b1 <= b2 && b1 < d)
                    {
                        D(p) = b1;
                        U(p) = std::max(u1,ix);
                        L(p) = std::min(l1,ix);
                    }
                    else
                    {
                        D(p) = b2;
                        U(p) = std::max(u2,ix);
                        L(p) = std::min(l2,ix);
                    }
                }
            }
        }
        else
        {
            uint8_t ix, u1, l1, b1, u2, l2, b2, d;

            for (i = height -1 ; i >0 ; i--)
            {
                for (j = width; j > 0; j--)
                {
                    point2d p = point2d(i,j);
                    ix = I(p);
                    u1 = U(point2d(i+1,j));
                    l1 = L(point2d(i+1,j));
                    b1 = std::max(u1,ix) - std::min(l1,ix);

                    u2 = U(point2d(i,j+1));
                    l2 = L(point2d(i,j+1));
                    b2 = std::max(u2,ix) - std::min(l2,ix);


                    d = D(p);



                    if (d <= b1 && d <= b2)
                        continue;
                    else if (b1 <= b2 && b1 < d)
                    {
                        D(p) = b1;
                        U(p) = std::max(u1,ix);
                        L(p) = std::min(l1,ix);
                    }
                    else
                    {
                        D(p) = b2;
                        U(p) = std::max(u2,ix);
                        L(p) = std::min(l2,ix);
                    }
                }
            }
        }
    }


    // Rasterscan


    double stop_s=clock();
    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

    io::imsave(D, "dkm.pgm");







}
