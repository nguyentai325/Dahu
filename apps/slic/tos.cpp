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
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>


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




void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb]   \n";
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
    image2d<V> Img;
    io::imread(fileName.c_str(), Img);
    box2d D = Img.domain();


    image2d<uint8_t> gray(D);

    mln_foreach(auto p, D)
    {
        gray(p) = 0.3 * Img(p)[0] + 0.59 *Img(p)[1] + 0.11*Img(p)[2];
    }


    // immerse

    box2d dom = gray.domain();
    dom.pmin = dom.pmin * 2;
    dom.pmax = dom.pmax * 2 - 1;
    image2d<std::pair<uint8_t,uint8_t>> out(dom);
    typedef point2d      P;
    mln_foreach(point2d p, gray.domain())
    {
        uint8_t a = gray.at(p),
          b = gray.at(p + P{0,1}),
          c = gray.at(p + P{1,0}),
          d = gray.at(p + P{1,1});

        uint8_t min1 = inf(a,b), min2 = inf(a,c);
        uint8_t max1 = sup(a,b), max2 = sup(a,c);
        uint8_t min3 = inf(d, inf(c, min1));
        uint8_t max3 = sup(d, sup(c, max1));

        point2d q = 2 * p;
        out.at(q) = std::make_pair (gray.at(p),gray.at(p));
        out.at(q + P{0,1}) = std::make_pair(min1, max1);
        out.at(q + P{1,0}) = std::make_pair(min2, max2);
        out.at(q + P{1,1}) = std::make_pair(min3, max3);
    }


    io::imsave(gray, "vkl.pgm");





}
