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
  std::cout << "Usage: " << argv[0] << " input[rgb] output[uint8] \n";
  std::exit(1);
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
}



int main(int argc, char** argv)
{
    if (argc < 3)
    usage(argv);
    const char* input_path = argv[1];
    const char* output_path = argv[2];

    std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_filter";
    std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/output_waterflow_mbd_gray/test_ima";
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

            double start_s=clock();


            image2d<V> Img;
            io::imread(fileName.c_str(), Img);
            box2d Dom = Img.domain();

            image2d<uint8>  ima_gray = rgb2gray(Img);

            image2d<uint8> I = addborder_gray(ima_gray);


            double stop_s=clock();
            std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

        //    # OUtput : MBD map D
        //    # Auxiliaries : U L M Q

        //    # State of the pixel
        //    # Set M(s) to flooded and initiate for s in S, other pixels are droughty


            // set of seeds S
            unsigned height = I.nrows();
            unsigned width = I.ncols();

            image2d<unsigned>  state(Dom);
            typedef std::pair<unsigned,unsigned> pair_t;
            image2d<pair_t> mm(Dom);
            image2d<uint8_t> dmap(Dom);
            //image2d<uint8_t> dmap1(Dom);

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
                        mm(p) = pair_t(I(p),I(p));
                    }
                    else
                    {
                        state(p) = 0;
                        dmap(p) = 255;
                        mm(p) = pair_t(I(p),I(p));
                    }
                }
            }

            stop_s=clock();
            std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;


            int dx[4] = {1 ,-1 , 0 , 0};
            int dy[4] = {0 , 0, 1, -1};

            // proceed the propagation of the pixel from border to the center of image until all of pixels is passed

            for (int lvl = 0; lvl < 256 ; lvl++)
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

                            if (state(r)==1 && dmap(r) > dmap(p))
                            {

                                mm(r) = mm(p);
                                if (I(r) < mm(r).first)
                                  mm(r).first = I(r);
                                if (I(r) > mm(r).second)
                                  mm(r).second = I(r);
                                if (dmap(r) > mm(r).second - mm(r).first)
                                {
                                    dmap(r) = mm(r).second - mm(r).first;
                                    Q[dmap(r)].push(r);
                                }
                            }

                            else if (state(r)==0)
                            {

                                mm(r) = mm(p);
                                if (I(r) < mm(r).first)
                                  mm(r).first = I(r);
                                if (I(r) > mm(r).second)
                                  mm(r).second = I(r);

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


            unsigned alpha = 30;
            unsigned beta = 200;


            mln_foreach(auto p, Dom)
            {

                if (dmap(p) <= alpha)
                    dmap(p) = 0;
                else if (dmap(p) >= beta)
                    dmap(p) = 255;
                else
                {
                    dmap(p) = uint8(float(dmap(p) - alpha)/ float(beta - alpha) *255);
                }

            }


            stop_s=clock();
            std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;


            fileName = outputDirectory + "/" + std::string(_dirent->d_name);
            io::imsave(dmap, fileName.c_str());

        }
    }

    closedir(directory);



}
