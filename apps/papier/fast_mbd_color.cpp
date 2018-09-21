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







void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] output[uint8]  \n";
  std::exit(1);
}



int main(int argc, char** argv)
{
    if (argc < 3)
    usage(argv);
    const char* input_path = argv[1];
    const char* output_path = argv[2];

    std::string inputDirectory = "/media/minh/DATA/Study/database/MSRA10K_Imgs_GT/MSRA10K_Imgs_GT/Imgs_resize";
    std::string outputDirectory = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/FastMBD/color";
    std::string outputDirectory_scalar = "/media/minh/DATA/Study/Results/compare_MBD_methods/MSRA/FastMBD/scalar";

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


            image2d<rgb8> I = addborder_color(Img);

            box2d Dom = I.domain();



            // number of passes K
            int K = 3;


            // set of seeds S
            unsigned height = I.nrows();
            unsigned width = I.ncols();

            std::cout << height<< std::endl;
            std::cout << width << std::endl;



            // MBD map



            image2d<rgb8> D(Dom);
            image2d<uint8> D_scalar(Dom);

            int i , j;

            for ( i= 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    if (i == 0 or i == height -1 or j == 0 or j == width -1 )
                        D(point2d(i,j)) = {0,0,0};
                    else
                        D(point2d(i,j)) = {255,255,255};
                }
            }



            // U,L

            image2d<rgb8> U(Dom);
            image2d<rgb8> L(Dom);


            mln_foreach(auto p , Dom)
            {
                U(p) = I(p);
                L(p) = I(p);
            }


            for (int k = 0; k < K; k++)
            {
                if (k % 2 == 0)
                {
                    rgb8 ix, u1, l1, b1, u2, l2, b2, d;

                    for (i = 1; i < height - 1 ; i++)
                    {
                        for (j = 1; j < width - 1; j++)
                        {
                            point2d p = point2d(i,j);
                            ix = I(p);
                            u1 = U(point2d(i-1,j));
                            l1 = L(point2d(i-1,j));
                            for (int t = 0 ; t < 3 ; t++)
                                b1[t] = std::max(u1[t],ix[t]) - std::min(l1[t],ix[t]);

                            u2 = U(point2d(i,j-1));
                            l2 = L(point2d(i,j-1));
                            for (int t = 0 ; t < 3 ; t++)
                                b2[t] = std::max(u2[t],ix[t]) - std::min(l2[t],ix[t]);


                            d = D(p);



                            if (d[0]+d[1]+d[2] <= b1[0]+b1[1]+b1[2] && d[0]+d[1]+d[2] <= b2[0] + b2[1] + b2[2])
                                continue;
                            else if (b1[0]+b1[1]+b1[2] <= b2[0] + b2[1] + b2[2] && b1[0]+b1[1]+b1[2] < d[0]+d[1]+d[2] )
                            {
                                D(p) = b1;
                                for (int t = 0 ; t < 3 ; t++)
                                {
                                    U(p)[t] = std::max(u1[t],ix[t]);
                                    L(p)[t] = std::min(l1[t],ix[t]);
                                }
                            }
                            else
                            {
                                D(p) = b2;
                                for (int t = 0; t < 3; t++)
                                {
                                    U(p)[t] = std::max(u2[t],ix[t]);
                                    L(p)[t] = std::min(l2[t],ix[t]);
                                }
                            }
                        }
                    }
                }
                else
                {
                    rgb8 ix, u1, l1, b1, u2, l2, b2, d;

                    for (i = height -1 ; i >0 ; i--)
                    {
                        for (j = width-1 ; j > 0; j--)
                        {
                            point2d p = point2d(i,j);
                            ix = I(p);
                            u1 = U(point2d(i+1,j));
                            l1 = L(point2d(i+1,j));
                            for(int t = 0 ; t <3; t++)
                                b1[t] = std::max(u1[t],ix[t]) - std::min(l1[t],ix[t]);

                            u2 = U(point2d(i,j+1));
                            l2 = L(point2d(i,j+1));
                            for(int t = 0 ; t <3; t++)
                                b2[t] = std::max(u2[t],ix[t]) - std::min(l2[t],ix[t]);


                            d = D(p);



                            if (d[0]+d[1]+d[2] <= b1[0]+b1[1]+b1[2] && d[0]+d[1]+d[2] <= b2[0] + b2[1] + b2[2])
                                continue;
                            else if (b1[0]+b1[1]+b1[2] <= b2[0] + b2[1] + b2[2] && b1[0]+b1[1]+b1[2] < d[0]+d[1]+d[2])
                            {
                                D(p) = b1;
                                for (int t = 0 ; t < 3 ; t++)
                                {
                                    U(p)[t] = std::max(u1[t],ix[t]);
                                    L(p)[t] = std::min(l1[t],ix[t]);
                                }
                            }
                            else
                            {
                                D(p) = b2;
                                for (int t = 0 ; t < 3 ; t++)
                                {
                                    U(p)[t] = std::max(u2[t],ix[t]);
                                    L(p)[t] = std::min(l2[t],ix[t]);
                                }
                            }
                        }
                    }
                }
            }


            unsigned alpha = 30;
            unsigned beta = 220;

            mln_foreach(auto p, U.domain())
            {

                D_scalar(p) = D(p)[0]/3 + D(p)[1]/3 + D(p)[2]/3;


            }

//            int max = 0;
//            mln_foreach(auto p, U.domain())
//            {
//                if (D_scalar(p) > max)
//                    max = D_scalar(p);
//            }

//            mln_foreach(auto p, U.domain())
//            {
//                D_scalar(p) = uint8(float(D_scalar(p))/ float(max) *255);
//            }

//            mln_foreach(auto p, U.domain())
//            {
//                if (D_scalar(p) <= alpha)
//                    D_scalar(p) = 0;
//                else if (D_scalar(p) >= beta)
//                    D_scalar(p) = 255;
//                else
//                {
//                    D_scalar(p) = uint8(float(D_scalar(p) - alpha)/ float(beta - alpha) *255);
//                }
//            }




            // Rasterscan


            double stop_s=clock();
            std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

            fileName = outputDirectory + "/" + std::string(_dirent->d_name);
            io::imsave(D, fileName.c_str());

            fileName = outputDirectory_scalar + "/" + std::string(_dirent->d_name);
            io::imsave(D_scalar, fileName.c_str());
            }

        }
    }

    closedir(directory);







}
