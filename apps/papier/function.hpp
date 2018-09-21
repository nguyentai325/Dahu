#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/dontcare.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/transpose_graph.hpp>
#include <boost/property_map/function_property_map.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/croutines.hpp>
#include <apps/g2/compute_g2.hpp>
#include <apps/g2/satmaxtree.hpp>

# include <mln/morpho/structural/dilate.hpp>


#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"


namespace mln
{


    ////////////////////////////// check border ////////////////////////////////////////
    /// \brief is_border check whether pixel belongs to the border
    /// \param p: pixel image
    /// \param height of image
    /// \param width of image
    /// \return true or false
    ///
    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 or  p[0] == height -1 or p[1] == 0  or p[1] == width -1 )
            return true;
        else
            return false;
    }

    ///////////////////////////// computed depth node on the tree /////////////////////////
    /// \brief depth_node_compute
    /// \param T: tree structure
    /// \return a vector of depth
    ///
    std::vector<unsigned>  depth_node_compute(tree_t T )
    {
        unsigned nodenumbers = T.nodes().size();
        std::vector<unsigned>  depth_node(nodenumbers);
        std::fill (depth_node.begin(),depth_node.end(),0);


        mln_foreach(auto x, T.nodes())
        {
            if (x.id() != 0)
                depth_node[x.id()] = depth_node[x.get_parent_id()] + 1;
        }
        return depth_node;

    }

    ////////////////////////// compute the Dahu distance between 2 points  ////////////////////////////
    /// \brief compute_dahu
    /// \param T: true structure
    /// \param depth_node: depth node vector on the tree
    /// \param node_value: image value of each node
    /// \param p1 : first pixel
    /// \param p2 : second pixel
    /// \return int: distance
    ///
    uint8 compute_dahu(tree_t T, std::vector<unsigned >depth_node, std::vector<rgb8> node_value, point2d p1, point2d p2)
    {

        auto& pmap  = T._get_data()->m_pmap;
        unsigned node_start = pmap(p1*2);
        unsigned node_dest = pmap(p2*2);



        unsigned node_cur = depth_node[node_start] >= depth_node[node_dest] ? node_start : node_dest;
        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest : node_start;

        std::vector<unsigned>  path;

        unsigned d_cur = depth_node[node_cur];
        unsigned d_obj = depth_node[node_obj];

        path.push_back(node_cur);

        //std::cout << "node object  "  << node_obj  << std::endl;



        rgb8 min, max;

        min = node_value[node_cur];
        max = node_value[node_cur];


        while (d_cur > d_obj)
        {
            d_cur -= 1;
            node_cur = T.get_node(node_cur).get_parent_id();
            path.push_back(node_cur);
            for (int i = 0; i < 3; i++)
            {
                if (node_value[node_cur][i] < min[i])
                  min[i] = node_value[node_cur][i];
                else if (node_value[node_cur][i] > max[i])
                  max[i] = node_value[node_cur][i];
            }
        }



        rgb8 d_dahu;
        uint8 d_dahu_scalar;

        while (node_cur != node_obj)
        {
            node_cur = T.get_node(node_cur).get_parent_id();
            node_obj = T.get_node(node_obj).get_parent_id();

            path.push_back(node_cur);
            path.push_back(node_obj);
            for (int i = 0; i < 3; i++)
            {
                if (node_value[node_cur][i] < min[i])
                  min[i] = node_value[node_cur][i];
                else if (node_value[node_cur][i] > max[i])
                  max[i] = node_value[node_cur][i];

                if (node_value[node_obj][i] < min[i])
                  min[i] = node_value[node_obj][i];
                else if (node_value[node_obj][i] > max[i])
                  max[i] = node_value[node_obj][i];
            }
        }


//        for (int i = 0; i < path.size(); i++)
//        {
//            std::cout << path[i]  << std::endl;
//            std::cout << int(node_value[path[i]][0])  << "  " <<  int(node_value[path[i]][1])  << "  " << int(node_value[path[i]][2])  << std::endl;
//        }
        d_dahu = max - min;
        d_dahu_scalar = d_dahu[0]/3 + d_dahu[1]/3 + d_dahu[2]/3;

        return d_dahu_scalar;
    }

    /////////////////////////// computed the water flow MBD distance between 2 points /////////////
    /// \brief compute_waterflow
    /// \param I : input image
    /// \param p_topleft : first pixel
    /// \param p_bottomright : second pixel
    /// \return int distance
    ///
    uint8 compute_waterflow(image2d<rgb8> I,  point2d p_topleft, point2d p_bottomright)
    {
        box2d Dom = I.domain();


        // set of seeds S
        unsigned height = I.nrows();
        unsigned width = I.ncols();

        image2d<unsigned>  state(Dom);
        typedef std::pair<rgb8,rgb8> pair_t;
        image2d<pair_t> mm(Dom);
        image2d<rgb8> dmap(Dom);

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
                if (i == p_topleft[0] and j == p_topleft[1] )
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

                if (p == p_bottomright)
                {
                    uint8 d_water_scalar = dmap(p_bottomright)[0]/3 + dmap(p_bottomright)[1]/3 + dmap(p_bottomright)[2]/3;
                    return d_water_scalar;
                }

                if (state(p) == 2)
                    continue;

                state(p) = 2;
                //dmap1(p) = lvl;

                for (int n1 = 0 ; n1 < 4 ; n1++)
                {
                    int x  = p[0] + dx[n1];
                    int y  = p[1] + dy[n1];

                    if (x >= 0 && x < height && y >= 0 and y < width)
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


    }


    ////////////////////////////// compute the MSTMBD between two points /////////////////////////////////////////////
    /// \brief compute_mstmbd
    /// \param T
    /// \param depth_node
    /// \param node_value
    /// \param p1
    /// \param p2
    /// \return
    ///

    uint8 compute_mstmbd(image2d<rgb8> f, image2d<point2d> parent, image2d<unsigned> depth_node_mst,  point2d p1, point2d p2)
    {

        point2d node_cur = depth_node_mst(p1) >= depth_node_mst(p2) ? p1 : p2;
        point2d node_obj = depth_node_mst(p1) >= depth_node_mst(p2) ? p2 : p1;

        unsigned d_cur = depth_node_mst(node_cur);
        unsigned d_obj = depth_node_mst(node_obj);


        rgb8 min, max;
        for (int i = 0 ; i < 3; i++)
        {
            min[i] = std::min(f(node_cur)[i],f(node_cur)[i]);
            max[i] = std::max(f(node_cur)[i],f(node_cur)[i]);
        }

        while (d_cur > d_obj)
        {
            d_cur -= 1;
            node_cur = parent(node_cur);
            for (int i = 0; i < 3; i++)
            {
                if (f(node_cur)[i] < min[i])
                  min[i] = f(node_cur)[i];
                else if (f(node_cur)[i] > max[i])
                  max[i] = f(node_cur)[i];
            }
        }

        rgb8 d_mstmbd;
        uint8 d_mstmbd_scalar;

        while (node_cur != node_obj)
        {
            node_cur = parent(node_cur);
            node_obj = parent(node_obj);

            for (int i = 0; i < 3; i++)
            {
                if (f(node_cur)[i] < min[i])
                  min[i] = f(node_cur)[i];
                else if (f(node_cur)[i] > max[i])
                  max[i] = f(node_cur)[i];

                if (f(node_obj)[i] < min[i])
                  min[i] = f(node_obj)[i];
                else if (f(node_obj)[i] > max[i])
                  max[i] = f(node_obj)[i];
                d_mstmbd[i] = max[i] - min[i];
            }
            d_mstmbd_scalar = d_mstmbd[0]/3 + d_mstmbd[1]/3 + d_mstmbd[2]/3;
        }
        return d_mstmbd_scalar;

    }


    //////////////////////////// compute the Dahu distance on the image //////////////////////
    /// \brief compute_dahu_image
    /// \param T : tree structure
    /// \param ima_compo : image value
    /// \return image_dahu
    /// "Saliency based detection of identity documents captured by smartphones"  Minh, Jonathan, Thierry Geraud
    ///
    image2d<rgb8> compute_dahu_image(tree_t T, image2d<rgb8> ima_compo)
    {


        unsigned nodenumbers = T.nodes().size();

        auto& pmap  = T._get_data()->m_pmap;

        box2d D = pmap.domain();
        unsigned height = D.pmax[0];
        unsigned width = D.pmax[1];

        image2d <rgb8>  distance_map(D);

        image2d <uint8_t>  distance_map_scalar(D);

        std::vector<rgb8>  min_node(nodenumbers);
        std::vector<rgb8>  max_node(nodenumbers);

        std::vector<rgb8>  min_temp_node(nodenumbers);
        std::vector<rgb8>  max_temp_node(nodenumbers);

        //std::vector<V>  dist(nodenumbers);
        std::vector<unsigned>  dist_gray_temp_node(nodenumbers);


        std::vector<int>  seed_node(nodenumbers);
        std::fill (seed_node.begin(),seed_node.end(),0);
        seed_node[0] = 1;
        std::vector<unsigned>  dejavu(nodenumbers);
        std::fill (dejavu.begin(),dejavu.end(),0);




        mln_foreach(auto p, D)
        {
            if (is_border(p,height,width))
            {
                dejavu[pmap(p)] = 1;
                seed_node[pmap(p)] = 1;
            }
        }






        mln_foreach(auto x, T.nodes())
        {


            std::vector<rgb8> listpixel1;

            mln_foreach (auto p, x.proper_pset())
            {
                listpixel1.push_back(ima_compo(ima_compo.point_at_index(p)));

            }


            std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end(),lexicographicalorder_less<value_t>());

            rgb8 medianvalue1;

            medianvalue1 = listpixel1[listpixel1.size()/2];
            //std::cout << "median value"  << medianvalue << "  " << medianvalue1  << std::endl;


            for (int i = 0; i< 3; ++i)
            {
                min_node[x.id()][i] = medianvalue1[i];
                max_node[x.id()][i] = medianvalue1[i];

                min_temp_node[x.id()][i] = medianvalue1[i];
                max_temp_node[x.id()][i] = medianvalue1[i];
                dist_gray_temp_node[x.id()] = 255 *3;


            }


            unsigned q = x.get_parent_id();

        }

        // chay up len


        std::cout << "chay up len  "  << std::endl;

        mln_reverse_foreach(auto x, T.nodes())
        {
            unsigned dis_temp = 0;
            rgb8 min_temp = min_node[x.get_parent_id()];
            rgb8 max_temp = max_node[x.get_parent_id()];

            if (dejavu[x.id()] == 1)
            {
                unsigned par_p = x.get_parent_id();

                if (seed_node[par_p] == 0)
                {

                    for (int i = 0; i< 3; ++i)
                    {
                        if (min_node[par_p][i] > min_temp_node[x.id()][i])
                            min_temp[i] = min_temp_node[x.id()][i];

                        if (max_node[par_p][i] < max_temp_node[x.id()][i])
                            max_temp[i] = max_temp_node[x.id()][i];

                        dis_temp = dis_temp + max_temp[i]  - min_temp[i];

                    }

                    if (dis_temp < dist_gray_temp_node[par_p])
                    {
                        dist_gray_temp_node[par_p] = dis_temp;
                        min_temp_node[par_p] = min_temp;
                        max_temp_node[par_p] = max_temp;
                    }
                    dejavu[par_p] =1;
                    std::cout << x.id() << "   " << par_p  << "   "<< dis_temp << std::endl;

                }
            }
        }

        // chay down xuong


        std::cout << "chay down xuong "  << std::endl;

        std::cout << "root " << T.get_root_id()  << std::endl;

        mln_foreach (auto x, T.nodes())
        {
            rgb8 min_temp = min_node[x.id()];
            rgb8 max_temp = max_node[x.id()];
            unsigned dis_temp = 0;

            if(seed_node[x.id()] == 0)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if (min_node[x.id()][i] > min_temp_node[x.get_parent_id()][i])
                        min_temp[i] = min_temp_node[x.get_parent_id()][i];
                    if (max_node[x.id()][i] < max_temp_node[x.get_parent_id()][i])
                        max_temp[i] = max_temp_node[x.get_parent_id()][i];
                    dis_temp = dis_temp + max_temp[i]  - min_temp[i];
                }

                if (dis_temp < dist_gray_temp_node[x.id()])
                {
                    dist_gray_temp_node[x.id()] = dis_temp;
                    max_temp_node[x.id()] = max_temp;
                    min_temp_node[x.id()] = min_temp;
                }
            }
        }


        // back propagation
        unsigned alpha = 20;
        unsigned beta = 160;

        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                distance_map(ima_compo.point_at_index(p)) = max_temp_node[x.id()] - min_temp_node[x.id()];
//                distance_map_scalar(ima_compo.point_at_index(p)) = distance_map(ima_compo.point_at_index(p))[0]/3 + distance_map(ima_compo.point_at_index(p))[1]/3 + distance_map(ima_compo.point_at_index(p))[2]/3;
//                if (distance_map_scalar(ima_compo.point_at_index(p)) <= alpha)
//                    distance_map_scalar(ima_compo.point_at_index(p)) = 0;
//                else if (distance_map_scalar(ima_compo.point_at_index(p)) >= beta)
//                    distance_map_scalar(ima_compo.point_at_index(p)) = 255;
//                else
//                {
//                    distance_map_scalar(ima_compo.point_at_index(p)) = uint8(float(distance_map_scalar(ima_compo.point_at_index(p)) - alpha)/ float(beta - alpha) *255);
//                }
            }
        }

        return distance_map;



    }



    /////////////////////////////////// compute the waterflow MBD on the image ///////////////////////////////////
    ///
    /// I input image
    /// return waterflow image
    /// "Water flow driven salient object detection at 180 fps" Huang, Zhang
    ///
    template <typename V>
    image2d<V> waterflow_mbd_distance(image2d<V> I)
    {
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

                    if (x >= 0 && x < height && y >= 0 and y < width)
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

        return dmap;

    }


    /////////////////////////////////// compute fast MBD distance image //////////////////////////
    /// I input image
    /// return : fast MBD
    /// "Minimum Barrier Distance object detection at 80 fps"  Zhang, sclaroff
    template <typename V>
    image2d<V> fast_mbd_distance(image2d<V> I)
    {

        box2d Dom = I.domain();

        // number of passes K
        int K = 3;


        // set of seeds S
        unsigned height = I.nrows();
        unsigned width = I.ncols();


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



        mln_foreach(auto p, U.domain())
        {
            D_scalar(p) = D(p)[0]/3 + D(p)[1]/3 + D(p)[2]/3;
        }

        return D;


    }


    /////////////////////////// construct a minimum spanning tree ////////////////////////////////////
    /// \brief primMST_color
    /// \param I input image
    /// \param S vector of pixels
    /// \return parent image
    ///
    image2d<point2d> primMST_color(image2d<rgb8> I, std::vector<point2d> &S, point2d p1)
    {
         image2d<point2d> parent(I.domain()); // Array to store constructed MST
         image2d<int> key(I.domain());   // Key values used to pick minimum weight edge in cut
         image2d<bool> mstSet(I.domain()); // To represent set of vertices not yet included in MST
         image2d<int> state(I.domain());

         unsigned height = I.nrows();
         unsigned width = I.ncols();


        mln_foreach(auto p, I.domain())
        {
             key(p) = INT_MAX;
             mstSet(p) = false;
             state(p) = 0;
         }

        typedef std::pair<int, point2d> iPair;

        //std::queue<point2d>  Q;
        std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > Q;

        point2d p = p1;
        key(p) = 0;
        parent(p) = p;
        Q.push(std::make_pair(0,p));


        int x[4] = {-1,1,0,0};
        int y[4] = {0,0,-1,1};



        while (!Q.empty())
        {
            point2d u = Q.top().second;
            Q.pop();

            if (mstSet(u) == true)
                    continue;

            mstSet(u) = true;
            S.push_back(u);

            for (int t = 0; t < 4; t++)
            {
                 int x_n = u[0] + x[t];
                 int y_n = u[1] + y[t];

                 if( (x_n >= 0 && x_n < height) && (y_n >= 0 && y_n < width) )
                 {
                     point2d v = point2d(x_n,y_n);
                     int temp = std::abs(int(I(v)[0])-int(I(u)[0])) + std::abs(int(I(v)[1])-int(I(u)[1])) + std::abs(int(I(v)[2])-int(I(u)[2]));

                     if (mstSet(v) == false and temp < key(v))
                     {
                         parent(v) = u;
                         key(v) = temp;
                         Q.push(std::make_pair(key(v),v));
                     }

                 }
            }
        }


        return parent;
    }


    /////////////////////////// construct a minimum spanning tree ////////////////////////////////////
    /// \brief primMST_color
    /// \param I input image
    /// \param S vector of pixels
    /// \return parent image
    ///
    image2d<point2d> primMST(image2d<uint8> I, std::vector<point2d> &S)
    {
         image2d<point2d> parent(I.domain()); // Array to store constructed MST
         image2d<int> key(I.domain());   // Key values used to pick minimum weight edge in cut
         image2d<bool> mstSet(I.domain()); // To represent set of vertices not yet included in MST
         image2d<int> state(I.domain());

         unsigned height = I.nrows();
         unsigned width = I.ncols();


        mln_foreach(auto p, I.domain())
        {
             key(p) = INT_MAX;
             mstSet(p) = false;
             state(p) = 0;
         }

        typedef std::pair<int, point2d> iPair;

        //std::queue<point2d>  Q;
        std::priority_queue< iPair, std::vector <iPair> , std::greater<iPair> > Q;

        point2d p = point2d(0,0);
        key(p) = 0;
        parent(p) = p;
        Q.push(std::make_pair(0,p));


        int x[4] = {-1,1,0,0};
        int y[4] = {0,0,-1,1};



        while (!Q.empty())
        {
            point2d u = Q.top().second;
            Q.pop();

            if (mstSet(u) == true)
                    continue;

            mstSet(u) = true;
            S.push_back(u);

            for (int t = 0; t < 4; t++)
            {
                 int x_n = u[0] + x[t];
                 int y_n = u[1] + y[t];

                 if( (x_n >= 0 && x_n < height) && (y_n >= 0 && y_n < width) )
                 {
                     point2d v = point2d(x_n,y_n);
                     int temp = std::abs(int(I(v))-int(I(u)));

                     if (mstSet(v) == false and temp < key(v))
                     {
                         parent(v) = u;
                         key(v) = temp;
                         Q.push(std::make_pair(key(v),v));
                     }

                 }
            }
        }


        return parent;
    }



    ////////////////////////////////  compute MBD distance image on the Minimum spanning tree //////////////////////////////
    /// I input image
    ///

    template <typename V>
    image2d<V> mbd_mst_distance(image2d<V> a)
    {
        int height = a.nrows();
        int width = a.ncols();

        //// Minimum spanning tree ///////////////////

        std::cout << "minimum spanning tree  "  << std::endl;
        std::vector<point2d> S;
        image2d<point2d> parent = primMST_color(a, S);
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
                    std::cout << S[i] << "   " << par_p  << "   "<< dis_temp << std::endl;
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


        unsigned alpha = 20;
        unsigned beta = 160;



        mln_foreach(auto p, a.domain())
        {
            for (int k = 0; k< 3; ++k)
            {
                distance_map(p)[k] = max_temp_node(p)[k] - min_temp_node(p)[k];
            }
            distance_map_scalar(p) = distance_map(p)[0]/3 + distance_map(p)[1]/3 + distance_map(p)[2]/3;
        }
        return distance_map;
    }




    ///////////////////////////// computed depth node on the tree /////////////////////////
    /// \brief depth_node_compute
    /// \param T: tree structure
    /// \return a vector of depth
    ///
    image2d<unsigned>  depth_node_mst_compute(image2d<point2d> parent, std::vector<point2d> S )
    {
        image2d<unsigned>  depth_node(parent.domain());
        depth_node(point2d(0,0)) = 0;

        for(int i = 0; i < S.size() ; i++)
        {
            if (S[i] != point2d(0,0))
                depth_node(S[i]) = depth_node(parent(S[i])) + 1;
        }
        return depth_node;

    }
    
    
    
    ////////////////////////////////////ROI //////////////////////////////////////////////////
    
        image2d<uint8_t> roi(tree_t T, image2d<rgb8> F, point2d p1, point2d p2)
    {

        unsigned nodenumbers = T.nodes().size();

        auto& nodes1 = T._get_data()->m_nodes;
        auto& S     = T._get_data()->m_S;
        auto& pmap  = T._get_data()->m_pmap;

        box2d D = pmap.domain();

        std::vector<unsigned>  depth_node(nodenumbers);
        std::fill (depth_node.begin(),depth_node.end(),0);


        mln_foreach(auto x, T.nodes())
        {
            if (x.id() != 0)
                depth_node[x.id()] = depth_node[x.get_parent_id()] + 1;
        }


        // find path and lca

        std::vector<unsigned> path;

        unsigned node_1 = pmap(p1);
        unsigned node_2 = pmap(p2);

        unsigned node_start = node_1;
        unsigned node_dest = node_2;

        unsigned node_cur = depth_node[node_start] >= depth_node[node_dest] ? node_start : node_dest;
        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest :node_start  ;

        unsigned d_cur = depth_node[node_cur];
        unsigned d_obj = depth_node[node_obj];


        while (d_cur > d_obj)
        {
            d_cur -= 1;
            path.push_back(node_cur);
            node_cur = T.get_node(node_cur).get_parent_id();
        }

        while (node_cur != node_obj)
        {
            path.push_back(node_cur);
            path.push_back(node_obj);
            node_cur = T.get_node(node_cur).get_parent_id();
            node_obj = T.get_node(node_obj).get_parent_id();
        }

        path.push_back(node_cur);

        unsigned lca = node_cur;


        // show tree path


        std::vector<uint8_t>  path_node(nodenumbers);
        std::fill (path_node.begin(),path_node.end(),0);

        for(int index = 0 ; index < path.size() ; ++index)
        {
            path_node[path[index]] = 255;
        }




        image2d<uint8_t> under_ima(D);
        extension::fill(under_ima, 0);

        std::vector<uint8_t>  under_node(nodenumbers);
        std::fill (under_node.begin(),under_node.end(),0);


        mln_foreach(auto x, T.nodes())
        {
            if (path_node[x.get_parent_id()] == 255 and path_node[x.id()] == 0)
            {
                under_node[x.id()] = 255;
            }
        }

        mln_foreach(auto x, T.nodes())
        {
            if (under_node[x.get_parent_id()] == 255 and under_node[x.id()] == 0)
            {
                under_node[x.id()] = 255;
            }
            mln_foreach (auto p, x.proper_pset())
            {
                under_ima(F.point_at_index(p)) = under_node[x.id()];
            }
        }


        image2d<uint8_t> res_ima(D);
        extension::fill(res_ima, 0);

        std::vector<uint8_t>  res_node(nodenumbers);
        std::fill (res_node.begin(),res_node.end(),0);
        res_node[lca] = 255;

        mln_foreach(auto x, T.nodes())
        {
            if (res_node[x.get_parent_id()] == 255 and res_node[x.id()] == 0)
            {
                res_node[x.id()] = 255;
            }
            mln_foreach (auto p, x.proper_pset())
            {
                res_ima(F.point_at_index(p)) = res_node[x.id()];
            }
        }


        rect2d r = make_rectangle2d(3, 3);
        res_ima = morpho::structural::dilate(res_ima, r);


        image2d <uint8_t>  result_ima(D);
        extension::fill(result_ima, 0);

        mln_foreach (auto p , D)
        {
            if (res_ima(p) != under_ima(p))
            {
                result_ima(p) =  255 ;
            }
        }


        return result_ima;
    }











}

#endif // FUNCTION_HPP
