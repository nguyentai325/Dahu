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

    std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_filter";
    std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/output_fastmbd_gray/test_ima";
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

            image2d<uint8>  ima_gray = rgb2gray(Img);

            image2d<uint8> I = addborder_gray(ima_gray);

            box2d Dom = I.domain();



            // number of passes K
            int K = 2;


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
                        for (j = width -1; j > 0; j--)
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


            unsigned alpha = 0;
            unsigned beta = 255;


//            mln_foreach(auto p, D.domain())
//            {

//                if (D(p) <= alpha)
//                    D(p) = 0;
//                else if (D(p) >= beta)
//                    D(p) = 255;
//                else
//                {
//                    D(p) = uint8(float(D(p) - alpha)/ float(beta - alpha) *255);
//                }

//            }
            int max = 0;
            mln_foreach(auto p, D.domain())
            {
                if (D(p) > max)
                    max = D(p);
            }

            mln_foreach(auto p, D.domain())
            {
                D(p) = uint8(float(D(p))/ float(max) *255);
            }



            // Rasterscan


            double stop_s=clock();
            std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

//            std::string fileout = output_path;
//            io::imsave(D, fileout.c_str());

            fileName = outputDirectory + "/" + std::string(_dirent->d_name);

            io::imsave(D, fileName.c_str());

        }
    }

    closedir(directory);





}
