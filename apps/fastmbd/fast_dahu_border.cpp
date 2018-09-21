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

#include "immerse.hpp"





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



    // load image
    image2d<V> Img;
    io::imread(fileName.c_str(), Img);
    box2d Dom = Img.domain();

    image2d<uint8_t> I(Dom);


    // convert rgb2gray image

    mln_foreach(auto p, Dom)
    {
        I(p) = 0.2989 * Img(p)[0] + 0.5870 * Img(p)[1] + 0.1140 * Img(p)[2];
    }


    // Preprocessing image (Boundary replace) =============

    image2d<uint8_t> pre_I(Dom);

    // compute average boundary
    int tong = 0;
    int dem = 1;
    mln_foreach(point2d p, I.domain())
    {

        if ( is_border(p , I))
        {
            tong = tong +  I(p);
            dem = dem +1;
        }
    }
    uint8_t average = tong/dem;

    std::cout << " average boundary  " <<int(average) << std::endl;

    // boundary weight
    float w = 0.5;
    mln_foreach(point2d p, I.domain())
    {

        if ( is_border(p , I))
            pre_I(p) = (1-w)*I(p) + w*average;
        else
            pre_I(p) = I(p);
    }



    // immerse image


    image2d< range<uint8_t> > U = immerse(pre_I);
    box2d d = U.domain();

    // set of seeds S
    unsigned height = U.nrows();
    unsigned width = U.ncols();

    image2d<unsigned>  state(d);
    typedef std::pair<unsigned,unsigned> pair_t;
    image2d<pair_t> mm(d);
    image2d<uint8_t> dmap(d);

    image2d<uint8_t> Ub(d);


    // priority queue

    std::vector<std::queue<point2d> > Q(256);



    // put seed on the border of the image
    // change the state of the pixel
    for (int i = 0; i < height ; i++)
    {
        for(int j = 0; j < width ; j++)
        {
            point2d p = point2d(i,j);
            if (i == 0 or i == height - 1 or j == 0 or j == width -1)
            {
                state(p) = 1;
                dmap(p) = 0;
                //dmap1(p) = 0;
                Q[dmap(p)].push(p);
                Ub(p) = U(p).lower;
                mm(p) = pair_t(Ub(p),Ub(p));
            }
            else
            {
                state(p) = 0;
                dmap(p) = 255;
                //mm(p) = pair_t(U(p),U(p));
            }
        }
    }


    int dx[4] = {1 ,-1 , 0 , 0};
    int dy[4] = {0 , 0, 1, -1};

    // proceed the propagation of the pixel from border to the center of image until all of pixels is passed

    for (int lvl = 0; lvl < 256 ; lvl++)
    {
        while (!Q[lvl].empty())
        {
            point2d p = Q[lvl].front();
            Q[lvl].pop();
            uint8_t l_cur = Ub(p);

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

                    //

                    uint8_t l_ ;
                    if (l_cur < U(r).lower)
                        l_ = U(r).lower;
                    else if (l_cur > U(r).upper)
                        l_ = U(r).upper;
                    else
                        l_ = l_cur;

                    Ub(r) = l_;

                    if (state(r)==1 && dmap(r) > dmap(p))
                    {

                        mm(r) = mm(p);

                        if (Ub(r) < mm(r).first)
                          mm(r).first = Ub(r);
                        if (Ub(r) > mm(r).second)
                          mm(r).second = Ub(r);
                        if (dmap(r) > mm(r).second - mm(r).first)
                        {
                            dmap(r) = mm(r).second - mm(r).first;
                            Q[dmap(r)].push(r);
                        }
                    }

                    else if (state(r)==0)
                    {

                        mm(r) = mm(p);

                        if (Ub(r) < mm(r).first)
                          mm(r).first = Ub(r);
                        if (Ub(r) > mm(r).second)
                          mm(r).second = Ub(r);

                        dmap(r) = mm(r).second - mm(r).first;
                        Q[dmap(r)].push(r);
                        state(r) =1;

                    }
                    else
                        continue;

                }

            }
        }
    }





    double stop_s=clock();
    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;




    io::imsave(dmap, "dmap_dahu_border.pgm");
    io::imsave(pre_I, "preprocessing.pgm");




}
