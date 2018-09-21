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
}



int main(int argc, char** argv)
{
    if (argc < 3)
    usage(argv);
    const char* input_path = argv[1];
    const char* output_path = argv[2];

    std::string inputDirectory = "/media/minh/DATA/Study/database/MSRA10K_Imgs_GT/MSRA10K_Imgs_GT/Imgs_resize";
    std::string outputDirectory = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/Waterflow/color";
    std::string outputDirectory_scalar = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/Waterflow/scalar";

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
            std::cout << fileName << std::endl;
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

            using namespace mln;
            typedef rgb8 V;

            double start_s=clock();


            image2d<V> Img;
            io::imread(fileName.c_str(), Img);

            image2d<V> I = addborder_color(Img);

            box2d Dom = I.domain();



        //    double stop_s=clock();
        //    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

        ////    # OUtput : MBD map D
        ////    # Auxiliaries : U L M Q

        ////    # State of the pixel
        ////    # Set M(s) to flooded and initiate for s in S, other pixels are droughty


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
                    if (i == 0 or i == height - 1 or j == 0 or j == width -1)
                    {
                        state(p) = 1;
                        dmap(p) = {0,0,0};
                        Q[dmap(p)[0]+dmap(p)[1]+dmap(p)[2]].push(p);
                        mm(p) = pair_t(I(p),I(p));
                    }
                    else
                    {
                        state(p) = 0;
                        dmap(p) = {255,255,255};
                        mm(p) = pair_t(I(p),I(p));
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

                            }
                            else
                                continue;

                        }

                    }
                }
            }


            unsigned alpha = 30;
            unsigned beta = 220;

            mln_foreach(auto p, Dom)
            {

                dmap_scalar(p) = dmap(p)[0]/3 + dmap(p)[1]/3 + dmap(p)[2]/3;
//                if (dmap_scalar(p) <= alpha)
//                    dmap_scalar(p) = 0;
//                else if (dmap_scalar(p) >= beta)
//                    dmap_scalar(p) = 255;
//                else
//                {
//                    dmap_scalar(p) = uint8(float(dmap_scalar(p) - alpha)/ float(beta - alpha) *255);
//                }

            }

//            int max = 0;
//            mln_foreach(auto p, Dom)
//            {
//                if (dmap_scalar(p) > max)
//                    max = dmap_scalar(p);
//            }

//            mln_foreach(auto p, Dom)
//            {
//                dmap_scalar(p) = uint8(float(dmap_scalar(p))/ float(max) *255);
//            }

//            mln_foreach(auto p, Dom)
//            {
//                if (dmap_scalar(p) <= alpha)
//                    dmap_scalar(p) = 0;
//                else if (dmap_scalar(p) >= beta)
//                    dmap_scalar(p) = 255;
//                else
//                {
//                    dmap_scalar(p) = uint8(float(dmap_scalar(p) - alpha)/ float(beta - alpha) *255);
//                }
//            }



            fileName = outputDirectory + "/" + std::string(_dirent->d_name);
            io::imsave(dmap, fileName.c_str());

            fileName = outputDirectory_scalar + "/" + std::string(_dirent->d_name);
            io::imsave(dmap_scalar, fileName.c_str());
            }

        }
    }

    closedir(directory);


}
