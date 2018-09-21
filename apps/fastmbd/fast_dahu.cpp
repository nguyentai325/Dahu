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






namespace mln
{

    uint8_t upper_level_next_to_lcur(std::vector<std::queue<point2d> >& q, uint8_t l_cur,
                                           bool& found)
    {
      uint8_t v = l_cur;
      for (;;)
        {
          if (! q[v].empty())
            {
              found = true;
              return v;
            }
          if (v == 255)
            break;
          v = v + 1;
        }
      found = false;
      return l_cur;
    }


    uint8_t lower_level_next_to_lcur(std::vector<std::queue<point2d> >& q, uint8_t l_cur,
                                           bool& found)
    {
      uint8_t v = l_cur;
      for (;;)
        {
          if (! q[v].empty())
            {
              found = true;
              return v;
            }
          if (v == 0)
            break;
          v = v - 1;
        }
      found = false;
      return l_cur;
    }


    uint8_t level_next_to_lcur(std::vector<std::queue<point2d> >& q, uint8_t& l_cur)
    {
      uint8_t l_;
      bool found;

      bool up = int(2. * std::rand() / (RAND_MAX + 1.));
      if (up)
        {
          l_ = upper_level_next_to_lcur(q, l_cur, found);
          if (found)
            return l_;
          else
            {
              l_ = lower_level_next_to_lcur(q, l_cur, found);
              if (! found)
                std::abort();
              return l_;
            }
        }
      else
        {
          l_ = lower_level_next_to_lcur(q, l_cur, found);
          if (found)
            return l_;
          else
            {
              l_ = upper_level_next_to_lcur(q, l_cur, found);
              if (! found)
                std::abort();
              return l_;
            }
        }
    }



    point2d priority_pop(std::vector<std::queue<point2d> >& q, uint8_t& l_cur, unsigned & n)
    // modify q, and sometimes l_cur
    {
        if (q[l_cur].empty())
        {
          uint8_t l_ = level_next_to_lcur(q, l_cur);  // such as q[l_] is not empty
          if (q[l_].empty())
            std::abort();
          l_cur = l_;
        }
        point2d r = q[l_cur].front();
        q[l_cur].pop();
        n = n -1;
        return r;

    }

    uint8_t
    priority_push(std::vector<std::queue<point2d> >& q, point2d p,
                  const image2d< range<uint8_t> >& U,
                  uint8_t l_cur, unsigned & n)
    // modify q
    {
      uint8_t
        lower = U(p).lower,
        upper = U(p).upper,
        l_;
      if (lower > l_cur)
        l_ = lower;
      else if (upper < l_cur)
        l_ = upper;
      else
        l_ = l_cur;
      q[l_].push(p);
      n = n +1;
      return l_;
    }

}



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

    mln_foreach(auto p, Dom)
    {
        I(p) = 0.2989 * Img(p)[0] + 0.5870 * Img(p)[1] + 0.1140 * Img(p)[2];
    }

    // immerse image


    image2d< range<uint8_t> > U = immerse(I);
    box2d d = U.domain();

    // propagate and compute dmap
    std::vector<std::queue<point2d> > q(256); // vector of queues

    image2d<bool> deja_vu(d);
    extension::fill(deja_vu, false);

    // distance map
    typedef std::pair<unsigned,unsigned> pair_t;
    image2d<pair_t> mm(d);

    image2d<uint8_t> dmap(d);

    point2d p_oo = point2d(0,0);

    uint8_t l_cur = U(p_oo).lower;

    q[l_cur].push(p_oo);
    deja_vu(p_oo) = true;


    mm(p_oo) = pair_t(l_cur, l_cur);
    dmap(p_oo) = 0;


    unsigned n = 1;


    int dx[4] = {1 ,-1 , 0 , 0};
    int dy[4] = {0 , 0, 1, -1};

    while(n != 0)
    {
        point2d p = priority_pop(q, l_cur, n);
        for (int n1 = 0 ; n1 < 4 ; n1++)
        {
            int x  = p[0] + dx[n1];
            int y  = p[1] + dy[n1];
            point2d nei = point2d(x,y);

            if (d.has(nei)  && deja_vu(nei) == false)
            {
                uint8_t l_ = priority_push(q, nei, U, l_cur,n);
                deja_vu(nei) = true;
                mm(nei) = mm(p);
                if (l_ < mm(nei).first)
                  mm(nei).first = l_;
                if (l_ > mm(nei).second)
                  mm(nei).second = l_;
                dmap(nei) = mm(nei).second - mm(nei).first;
            }

        }

    }




    double stop_s=clock();
    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;



    io::imsave(dmap, "dkmm.pgm");





}
