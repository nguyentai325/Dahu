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
#include <mln/morpho/structural/dilate.hpp>

#include "immerse.hpp"





#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"
#include "function.hpp"



// Compute the depth attribute of each graph node
boost::vector_property_map<unsigned>
compute_graph_depth(const MyGraph& g)
{
  mln_entering("Compute graph depth");

  boost::vector_property_map<unsigned> depth(boost::num_vertices(g));

  auto one = [](mln::dontcare_t) -> int{ return 1; };
  auto w = boost::make_function_property_map<MyGraph::edge_descriptor, int, decltype(one)>(one);

  MyGraph::vertex_descriptor root = boost::vertex(0, g);
  depth[root] = 0;

  MyGraph gT;
  boost::transpose_graph(g, gT);
  boost::dag_shortest_paths(gT, root, boost::weight_map(w)
                            .distance_map(depth)
                            .distance_compare(std::greater<int> ())
                            .distance_inf(-1)
                            .distance_zero(0)
                            );
  mln_exiting();
  return depth;
}

// Compute the per-pixel attribute and reconstruct
template <class ValueMap>
void
write_vmap_to_image(const tree_t* t, const tlink_t* tlink,
                    const ValueMap& vmap, mln::image2d<mln::uint16>& out)
{
  mln_foreach(auto px, out.pixels())
  {
    unsigned w = 0;
    for (int i = 0; i < NTREE; ++i)
      {
        tree_t::node_type tnode = t[i].get_node_at(px.index());
        MyGraph::vertex_descriptor gnode = tlink[i][tnode];
        w = std::max(w, vmap[gnode]);
      }
    px.val() = w;
  }
}

/// \brief Remove non-2F from the tree
template <class P>
mln::morpho::component_tree<P, mln::image2d<P> >
tree_keep_2F(const mln::morpho::component_tree<P, mln::image2d<P> >& tree)
{
  using namespace mln;
  morpho::component_tree<P, image2d<P> > out;

  auto newdata = out._get_data();
  auto olddata = tree._get_data();

  // 1. Copy the point2node map
  box2d olddom = olddata->m_pmap.domain();
  box2d dom;
  dom.pmin = olddom.pmin / 2;
  dom.pmax = (olddom.pmax + 1) / 2;
  newdata->m_pmap.resize(dom);
  copy(olddata->m_pmap | sbox2d{olddom.pmin, olddom.pmax, {2,2}},
       newdata->m_pmap);

  // 2. Copy the node
  newdata->m_nodes = olddata->m_nodes;

  // 3. Copy the point set and update node first point index/
  newdata->m_S.resize(dom.size());
  unsigned j = 0;
  for (unsigned i = 0; i < olddata->m_S.size(); ++i)
    {
      P p = olddata->m_S[i];
      point2d pt = olddata->m_pmap.point_at_index(p);
      if (K1::is_face_2(pt))
        {
          newdata->m_S[j] = newdata->m_pmap.index_of_point(pt/2);
          auto node = tree.get_node_at(p);
          if (node.get_first_point_id() == i)
            newdata->m_nodes[node.id()].m_point_index = j;
          ++j;
        }
    }
  // 4. Do not forget the sentinel
  newdata->m_nodes[out.npos()].m_point_index = j;

  return out.get_subtree(tree.get_root_id());
}



namespace mln
{

//    bool is_border (point2d p, unsigned height , unsigned width)
//    {
//        if (p[0] == 0 or  p[0] == height -1 or p[1] == 0  or p[1] == width -1 )
//            return true;
//        else
//            return false;
//    }

//    bool is_border (point2d p, unsigned height , unsigned width)
//    {
//        if (p[0] != 0  and p[0]  != 1  and p[0] != height -1 and p[0] != height -2 and p[1] != 0  and p[1]  != 1 and p[1]  != width - 1 and p[1] != width -2)
//        {
//            if (p[0] == 2 or p[0] == height - 3 or p[1] == width -3 or p[1] == 2  )
//                return true;
//            else
//                return false;
//        }
//    }

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

    template <typename V>
    image2d <V> dahu_distance_2_points(tree_t T, image2d<V> ima_compo, point2d p_topleft, point2d p_bottomright, V & dahu_distance )
    {
        // Initiation
        unsigned height = ima_compo.nrows();
        unsigned width = ima_compo.ncols();

        unsigned nodenumbers = T.nodes().size();

        auto& nodes1 = T._get_data()->m_nodes;
        auto& S     = T._get_data()->m_S;
        auto& pmap  = T._get_data()->m_pmap;

        box2d D = pmap.domain();

        image2d <V>  distance_map(D);

        image2d <uint8_t>  distance_map_scalar(D);

        std::vector<V>  min_node(nodenumbers);
        std::vector<V>  max_node(nodenumbers);

        std::vector<V>  min_temp_node(nodenumbers);
        std::vector<V>  max_temp_node(nodenumbers);

        std::vector<V>  dist(nodenumbers);
        std::vector<unsigned>  dist_gray_temp_node(nodenumbers);


        std::vector< std::vector<unsigned> > children(nodenumbers);

        std::vector<int>  seed_node(nodenumbers);
        std::fill (seed_node.begin(),seed_node.end(),0);
        seed_node[0] = 1;
        std::vector<unsigned>  dejavu(nodenumbers);
        std::fill (dejavu.begin(),dejavu.end(),0);







        mln_foreach(auto p, D)
        {
//            if (is_border(p,height,width))
//            {
//                dejavu[pmap(p)] = 1;
//                seed_node[pmap(p)] = 1;
//            }
            if(p == p_topleft)
            {
                dejavu[pmap(p)] = 1;
                seed_node[pmap(p)] = 1;

            }
        }




        std::vector<V>  color_node(nodenumbers);


        mln_foreach(auto x, T.nodes())
        {


            std::vector<V> listpixel1;

            mln_foreach (auto p, x.proper_pset())
            {
                listpixel1.push_back(ima_compo(ima_compo.point_at_index(p)));

            }


            std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end(),lexicographicalorder_less<value_t>());

            V medianvalue1;

            medianvalue1 = listpixel1[listpixel1.size()/2];

            for (int i = 0; i< 3; ++i)
            {
                min_node[x.id()][i] = medianvalue1[i];
                max_node[x.id()][i] = medianvalue1[i];

                min_temp_node[x.id()][i] = medianvalue1[i];
                max_temp_node[x.id()][i] = medianvalue1[i];
                dist_gray_temp_node[x.id()] = 255 *3;


            }


            color_node[x.id()] = medianvalue1;

            unsigned q = x.get_parent_id();
            children[q].push_back(x.id());

        }


        // compute depth tree

        std::vector<unsigned>  depth_node(nodenumbers);
        std::fill (depth_node.begin(),depth_node.end(),0);


        mln_foreach(auto x, T.nodes())
        {
            if (x.id() != 0)
                depth_node[x.id()] = depth_node[x.get_parent_id()] + 1;
        }


        // find path and lca



        std::vector<unsigned> path;
        V Dahu_tree;
        V min_node_tree;
        V max_node_tree;

        unsigned node_topleft = pmap(p_topleft);
        unsigned node_bottomright = pmap(p_bottomright);


        unsigned node_start = node_topleft;
        unsigned node_dest = node_bottomright;


        unsigned node_cur = depth_node[node_start] >= depth_node[node_dest] ? node_start : node_dest;
        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest :node_start  ;



        unsigned d_cur = depth_node[node_cur];
        unsigned d_obj = depth_node[node_obj];


        //std::cout << d_cur   << "  "  << d_obj  << std::endl;

        min_node_tree = color_node[node_cur];
        max_node_tree = color_node[node_cur];

        //std::cout << "ola  "  << std::endl;
        while (d_cur > d_obj)
        {
            d_cur -= 1;
            path.push_back(node_cur);

            //std::cout << int(color_node[node_cur][0]) << "  " << int(color_node[node_cur][1]) << "  " << int(color_node[node_cur][2]) << std::endl;

            for (int k = 0 ; k < 3 ; k++)
            {
                if (color_node[node_cur][k] <  min_node_tree[k])
                    min_node_tree[k] = color_node[node_cur][k];
                if (color_node[node_cur][k] >  max_node_tree[k])
                    max_node_tree[k] = color_node[node_cur][k];
            }

            node_cur = T.get_node(node_cur).get_parent_id();
        }



        while (node_cur != node_obj)
        {
            path.push_back(node_cur);
            //std::cout << int(color_node[node_cur][0]) << "  " << int(color_node[node_cur][1]) << "  " << int(color_node[node_cur][2]) << std::endl;

            for (int k = 0 ; k < 3 ; k++)
            {
                if (color_node[node_cur][k] <  min_node_tree[k])
                    min_node_tree[k] = color_node[node_cur][k];
                if (color_node[node_cur][k] >  max_node_tree[k])
                    max_node_tree[k] = color_node[node_cur][k];
            }
            path.push_back(node_obj);
            //std::cout << int(color_node[node_obj][0]) << "  " << int(color_node[node_obj][1]) << "  " << int(color_node[node_obj][2]) << std::endl;

            for (int k = 0 ; k < 3 ; k++)
            {
                if (color_node[node_obj][k] <  min_node_tree[k])
                    min_node_tree[k] = color_node[node_obj][k];
                if (color_node[node_obj][k] >  max_node_tree[k])
                    max_node_tree[k] = color_node[node_obj][k];
            }
            node_cur = T.get_node(node_cur).get_parent_id();
            node_obj = T.get_node(node_obj).get_parent_id();
        }


        path.push_back(node_cur);

        for (int k = 0 ; k < 3 ; k++)
        {
            if (color_node[node_cur][k] <  min_node_tree[k])
                min_node_tree[k] = color_node[node_cur][k];
            if (color_node[node_cur][k] >  max_node_tree[k])
                max_node_tree[k] = color_node[node_cur][k];
            Dahu_tree[k] = max_node_tree[k] - min_node_tree[k];
        }
        unsigned lca = node_cur;
        dahu_distance = Dahu_tree;



        // chay up len


        std::cout << "chay up len  "  << std::endl;

        mln_reverse_foreach(auto x, T.nodes())
        {
            unsigned dis_temp = 0;
            V min_temp = min_node[x.get_parent_id()];
            V max_temp = max_node[x.get_parent_id()];

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
            V min_temp = min_node[x.id()];
            V max_temp = max_node[x.id()];
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
        unsigned alpha = 30;
        unsigned beta = 200;

        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                distance_map(ima_compo.point_at_index(p)) = max_temp_node[x.id()] - min_temp_node[x.id()];
                distance_map_scalar(ima_compo.point_at_index(p)) = distance_map(ima_compo.point_at_index(p))[0]/3 + distance_map(ima_compo.point_at_index(p))[1]/3 + distance_map(ima_compo.point_at_index(p))[2]/3;
                if (distance_map_scalar(ima_compo.point_at_index(p)) <= alpha)
                    distance_map_scalar(ima_compo.point_at_index(p)) = 0;
                else if (distance_map_scalar(ima_compo.point_at_index(p)) >= beta)
                    distance_map_scalar(ima_compo.point_at_index(p)) = 255;
                else
                {
                    distance_map_scalar(ima_compo.point_at_index(p)) = uint8(float(distance_map_scalar(ima_compo.point_at_index(p)) - alpha)/ float(beta - alpha) *255);
                }
            }
        }


        return distance_map;

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


    std::vector<point2d> tracePath(image2d<point2d> parent, point2d p_bottomright)
    {
        printf ("\nThe Path is ");
//        int row = dest.first;
//        int col = dest.second;
        point2d p = p_bottomright;

        std::vector<point2d> Path;

        while (!(parent(p) == p ))
        {
            Path.push_back(p);
            point2d temp_p = parent(p);
            p = temp_p;
        }

        Path.push_back(p);
//        while (!Path.empty())
//        {
//            point2d p = Path.top();
//            Path.pop();
//            printf("-> (%d,%d) ",p[0],p[1]);
//        }
        return Path;

    }

    std::vector<point2d> tracePath_mst(image2d<point2d> parent, point2d p_bottomright)
    {
        printf ("\nThe Path is ");
//        int row = dest.first;
//        int col = dest.second;
        point2d p = p_bottomright;

        std::vector<point2d> Path;

        while (!(parent(p) == p ))
        {
            Path.push_back(p);
            point2d temp_p = parent(p);
            p = temp_p;
        }

        Path.push_back(p);
//        while (!Path.empty())
//        {
//            point2d p = Path.top();
//            Path.pop();
//            printf("-> (%d,%d) ",p[0],p[1]);
//        }
        return Path;

    }



    template <typename V>
    image2d<V> waterflow_mbd_distance_2_points(image2d<V> I, point2d p_topleft, point2d p_bottomright, V & waterflow_distance, image2d<point2d> & parent )
    {
        box2d Dom = I.domain();

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
                if (p == p_topleft)
                {
                    state(p) = 1;
                    dmap(p) = 0;
                    //dmap1(p) = 0;
                    Q[dmap(p)].push(p);
                    mm(p) = pair_t(I(p),I(p));
                    parent(p) = p;
                }
                else
                {
                    state(p) = 0;
                    dmap(p) = 255;
                    mm(p) = pair_t(I(p),I(p));
                    parent(p) = point2d(-1,-1);
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
                                parent(r) = p;
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
                            parent(r) =p;

                        }
                        else
                            continue;

                    }

                }
            }
        }
        waterflow_distance  = dmap(p_bottomright);

        return dmap;

    }

    template <typename V>
    image2d<V> fast_dahu_2_points(image2d<V> I, point2d p_topleft, point2d p_bottomright, V & Dahuflow_distance, image2d<point2d> & parent )
    {
        // immerse image


        p_topleft = p_topleft *2 ;
        p_bottomright = p_bottomright *2;
        image2d< range<uint8_t> > U = immerse(I);
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
                if (p == p_topleft)
                {
                    state(p) = 1;
                    dmap(p) = 0;
                    //dmap1(p) = 0;
                    Q[dmap(p)].push(p);
                    Ub(p) = U(p).lower;
                    mm(p) = pair_t(Ub(p),Ub(p));
                    parent(p) = p;
                }
                else
                {
                    state(p) = 0;
                    dmap(p) = 255;
                    //mm(p) = pair_t(U(p),U(p));
                    parent(p) = point2d(-1,-1);
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
                                parent(r) = p;
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
                            parent(r) = p;

                        }
                        else
                            continue;

                    }

                }
            }
        }
        Dahuflow_distance  = dmap(p_bottomright);

        return dmap;


    }

    std::vector<point2d> tracePath_dahu(image2d<point2d> parent, point2d p_bottomright)
    {
        printf ("\nThe Path is ");
//        int row = dest.first;
//        int col = dest.second;
        point2d p = p_bottomright*2;

        std::vector<point2d> Path;

        while (!(parent(p) == p ))
        {
            Path.push_back(p);
            point2d temp_p = parent(p);
            p = temp_p;
        }

        Path.push_back(p);
//        while (!Path.empty())
//        {
//            point2d p = Path.top();
//            Path.pop();
//            printf("-> (%d,%d) ",p[0],p[1]);
//        }
        return Path;

    }



    image2d<uint8_t> roi(tree_t T, image2d<uint8> F, point2d p1, point2d p2)
    {

        p1 = p1 *2;
        p2 = p2 *2;

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


//        rect2d r = make_rectangle2d(3, 3);
//        res_ima = morpho::structural::dilate(res_ima, r);


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

    image2d<uint8_t> roi_alt(tree_t T, image2d<uint8> F, point2d p1, point2d p2)
    {

        p1 = p1 *2 ;
        p2 = p2 *2;
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
        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest : node_start;

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

        std::vector<unsigned> path_1;
        std::vector<unsigned> path_2;

        while (node_start != 0)
        {
            path_1.push_back(node_start);
            node_start = T.get_node(node_start).get_parent_id();
        }

        while (node_dest != 0)
        {
            path_2.push_back(node_dest);
            node_dest = T.get_node(node_dest).get_parent_id();
        }

        std::cout << "path  "  << std::endl;

        for(int index = 0 ; index < path.size() ; ++index)
        {
            std::cout << path[index]  << "   ";
        }

        std::cout << "  " << std::endl;


        std::cout << "path_1  "  << std::endl;

        for(int index = 0 ; index < path_1.size() ; ++index)
        {
            std::cout << path_1[index]  << "   ";
        }

        std::cout << "  " << std::endl;

        std::cout << "path_2  "  << std::endl;

        for(int index = 0 ; index < path_2.size() ; ++index)
        {
            std::cout << path_2[index]  << "   ";
        }

        std::cout << "  " << std::endl;


        // show tree path


        std::vector<uint8_t>  path_node(nodenumbers);
        std::fill (path_node.begin(),path_node.end(),0);

        for(int index = 0 ; index < path.size() ; ++index)
        {
            path_node[path[index]] = 255;
        }




        image2d<uint8_t> result_ima(D);
        extension::fill(result_ima, 0);



        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                result_ima(F.point_at_index(p)) = path_node[x.id()];
            }
        }

        return result_ima;
    }



    image2d<uint8_t> roi_alt_new(tree_t T, image2d<uint8> F, point2d p1, point2d p2)
    {

        p1 = p1 *2;
        p2 = p2 *2;
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
        unsigned node_obj = depth_node[node_start] >= depth_node[node_dest] ? node_dest : node_start;

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




        image2d<uint8_t> outpath_ima(D);
        extension::fill(outpath_ima, 0);

        image2d<uint8_t> result_ima(D);
        extension::fill(result_ima, 0);



        mln_foreach(auto x, T.nodes())
        {
            mln_foreach (auto p, x.proper_pset())
            {
                outpath_ima(F.point_at_index(p)) = 255 - path_node[x.id()];
            }
        }

        mln_foreach (auto p , D)
        {
            result_ima(p) = 255 - outpath_ima(p);
        }

        return result_ima;
    }




}



void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] α₀ α₁ λ output[rgb]\n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}


int main(int argc, char** argv)
{
  if (argc < 5)
    usage(argv);



  const char* input_path = argv[1];
  int a0 = std::atoi(argv[2]);
  int a1 = std::atoi(argv[3]);
  int lambda = std::atoi(argv[4]);
  const char* output_path = argv[5];


  tbb::task_scheduler_init init;


  std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_image";
  std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/shortest_path/test_ima";
  //std::string outputDirectory_scalar = "/home/minh/Documents/code/code_edwin/build/apps/papier/output_dahu_color/scalar";

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


  //std::string fileName = input_path;
  std::string fileName = inputDirectory + "/" +std::string(_dirent->d_name);


  //std::string fileName = input_path;

  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<rgb8> ima;
  io::imread(fileName.c_str(), ima);


  image2d<rgb8> pathimage(ima.domain());
  image2d<rgb8> pathimage_mst(ima.domain());
  image2d<rgb8> pathimage_dahu = interpolate_k1(ima);

  mln_foreach( auto p, ima.domain())
  {
      pathimage(p) = ima(p);
      pathimage_mst(p) = ima(p);
  }


  //image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
  image2d<uint8>  ima_gray = rgb2gray(ima);
  //io::imsave(ima_gray, "image_ori.pgm");
  image2d<uint8> f = ima_gray;
  image2d<uint8> F = interpolate_k1(f);
  unsigned height_F = F.nrows();
  unsigned width_F = F.ncols();

//  //////////////////////////////////   Tree of shapes ///////////////////////////////////////////////////////////////

    tree_t T;
    T = morpho::cToS(f, c4);
    if (a0 > 0)
    {
    grain_filter_inplace(T, a0);
    T.shrink_to_fit();
    }
    T._reorder_pset();

    auto& parent_pixel  = T._get_data()->m_parent_pixel;

    std::cout << parent_pixel.domain() << std::endl;





//    // 2 points

////    point2d p_topleft = point2d(130,130);
////    point2d p_bottomright = point2d(height_F-79,width_F-79);
////    point2d p_center = point2d(272, 364);

    //point2d p_topleft = point2d(1,1);
    point2d p_topleft = point2d(85,77);
    //point2d p_bottomright = point2d(height_F-35,width_F-35);
    //point2d p_center = point2d(100, 200);
    point2d p_center = point2d(60, 255);

    for (int i = 0; i < pathimage_dahu.nrows(); i++)
    {
        for(int j = 0; j < pathimage_dahu.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0]*2 - 4 and i < p_topleft[0]*2 + 4 and j > p_topleft[1]*2 - 4 and j < p_topleft[1]*2 + 4)  or (i > p_center[0]*2-4 and i < p_center[0]*2+4 and j > p_center[1]*2-4 and j < p_center[1]*2+4) )
            {
                pathimage_dahu(point2d(i,j)) = {255,0,0};
            }

        }
    }

    std::vector<point2d>  shortestPath_dahu1 = tracePath_dahu (parent_pixel, p_center);

    //image2d<uint8_t> pathimage_dahu(F.domain());

    for (int j = 0 ; j < shortestPath_dahu1.size() ; j++)
    {
        point2d p = shortestPath_dahu1[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage_dahu(p) = {0,0,255};
    }
    fileName = outputDirectory + "/"+ "dahu_" +std::string(_dirent->d_name);

    io::imsave(pathimage_dahu, fileName.c_str());


////    V dahu_distance;
////    image2d<V> dahu_map =  dahu_distance_2_points( T,  ima_compo,  p_topleft,  p_bottomright, dahu_distance );
////    io::imsave(dahu_map, "dmap1.png");

////    std::cout << p_bottomright << std::endl;
////    std::cout << int(dahu_map(p_bottomright)[0])  << "  " << int(dahu_map(p_bottomright)[1]) << "  "  <<int(dahu_map(p_bottomright)[2])  << std::endl;
////    std::cout << int(dahu_map(p_center)[0])  << "  " << int(dahu_map(p_center)[1]) << "  "  <<int(dahu_map(p_center)[2])  << std::endl;


//    ///////////////////////////////////////  Waterflow MBD  ////////////////////////////////////////////////////////////


    image2d<uint8> I = ima_gray;

    uint8 waterflow_distance;
    box2d Dom = I.domain();
    image2d<point2d > parent(Dom);
    image2d<uint8> waterflow_map = waterflow_mbd_distance_2_points( I,  p_topleft,  p_center, waterflow_distance, parent );
    //io::imsave(waterflow_map, "dmap2.png");
    std::vector<point2d>  shortestPath = tracePath (parent, p_center);



    for (int i = 0; i < ima.nrows(); i++)
    {
        for(int j = 0; j < ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_center[0]-4 and i < p_center[0]+4 and j > p_center[1]-4 and j < p_center[1]+4) )
            {
                pathimage(point2d(i,j)) = {255,0,0};
            }

        }
    }


    for (int j = 0 ; j < shortestPath.size() ; j++)
    {
        point2d p = shortestPath[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage(p)  = {0,0,255};
    }

//    // interpolation

    box2d D = Dom;
//    D.pmin = D.pmin * 2;
//    D.pmax = D.pmax * 2 - 1;

//    image2d<uint8_t> path_waterflow(D);
//    typedef point2d      P;

//    mln_foreach(point2d p, pathimage.domain())
//      {

//        uint8_t a = pathimage(p),
//          b = pathimage(p + P{0,1}),
//          c = pathimage(p + P{1,0}),
//          d = pathimage(p + P{1,1});


//        point2d q = 2 * p;
//        path_waterflow(q) = pathimage(p);

//        if (a == 255 and b == 255)
//            path_waterflow(q + P{0,1}) = 255;
//        else
//            path_waterflow(q + P{0,1}) = 0;

//        if (a == 255 and c == 255)
//            path_waterflow(q + P{1,0}) = 255;
//        else
//            path_waterflow(q + P{1,0}) = 0;

//        if (a == 255 and d == 255)
//            path_waterflow(q + P{1,1}) = 255;
//        else
//            path_waterflow(q + P{1,1}) = 0;

//      }


    fileName = outputDirectory + "/"+ "water_" +std::string(_dirent->d_name);

    io::imsave(pathimage, fileName.c_str());


    //    ///////////////////////////////////////  Minimum spanning tree MBD  ////////////////////////////////////////////////////////////

    std::vector<point2d> S1;
    image2d<point2d> parent_mst = primMST(f, S1);


    image2d<unsigned> depth_node_mst = depth_node_mst_compute(parent_mst, S1 );


    std::vector<point2d>  shortestPath_mst = tracePath_mst (parent_mst, p_center);



    for (int i = 0; i < ima.nrows(); i++)
    {
        for(int j = 0; j < ima.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_center[0]-4 and i < p_center[0]+4 and j > p_center[1]-4 and j < p_center[1]+4) )
            {
                pathimage_mst(point2d(i,j)) = {255,0,0};
            }

        }
    }

    for (int j = 0 ; j < shortestPath_mst.size() ; j++)
    {
        point2d p = shortestPath_mst[j];
        //printf("-> (%d,%d) ",p[0],p[1]);
        pathimage_mst(p) = {0,0,255};
    }
    fileName = outputDirectory + "/"+ "mst_" +std::string(_dirent->d_name);

    io::imsave(pathimage_mst, fileName.c_str());


//    //////////////////////////////////////////////////  Fast Dahu  //////////////////////////////////////


//    uint8 dahuflow_distance;
//    image2d<point2d > parent_dahu(D);
//    image2d<uint8> Dahuflow_map = fast_dahu_2_points( I,  p_topleft,  p_center, dahuflow_distance, parent_dahu );
//    io::imsave(Dahuflow_map, "dmap_dahu.png");

//    std::vector<point2d>  shortestPath_dahu = tracePath_dahu(parent_dahu, p_center);

//    image2d<uint8_t> pathimage_fastdahu(D);

//    for (int j = 0 ; j < shortestPath_dahu.size() ; j++)
//    {
//        point2d p = shortestPath_dahu[j];
//        //printf("-> (%d,%d) ",p[0],p[1]);
//        pathimage_fastdahu(p) = 255;
//    }
//    fileName = outputDirectory + "/"+ "fastdahu_" +std::string(_dirent->d_name);


//    io::imsave(pathimage_fastdahu, fileName.c_str());


////    box2d D1;
////    D1.pmin = point2d(0,0);
////    D1.pmax = point2d(9,9);
////    std::cout << D1 << std::endl;

////    image2d<uint8_t > test_image(D1);
////    for (int i = 0; i < I.nrows(); i++)
////    {
////        for(int j = 0; j < I.ncols() ; j++)
////        {
////            point2d p = point2d(i,j);
////            if (i >= 88/2 and i <= 106/2 and j >= 198/2 and j <= 216/2)
////            {
////                std::cout << i << "  "  << j << std::endl;
////                test_image(point2d(i-44,j-99)) = I(point2d(i,j));
////                //std::cout << int(test_image(point2d(i,j)))  << std::endl;
////            }
////        }
////    }
////    io::imsave(test_image, "image_test.pgm");

//    for (int i = 0; i < I.nrows(); i++)
//    {
//        for(int j = 0; j < I.ncols() ; j++)
//        {
//            point2d p = point2d(i,j);
//            if ((i > p_topleft[0] - 4 and i < p_topleft[0] + 4 and j > p_topleft[1] - 4 and j < p_topleft[1] + 4)  or (i > p_bottomright[0]-4 and i < p_bottomright[0]+4 and j > p_bottomright[1]-4 and j < p_bottomright[1]+4) )
//            {
//                I(point2d(i,j)) = 255;
//            }

//            if ((i > p_center[0]-4 and i < p_center[0]+4 and j > p_center[1]-4 and j < p_center[1]+4)  )
//            {
//                I(point2d(i,j)) = 255;
//            }
//        }
//    }
//    io::imsave(I, "image_ori.pgm");



//    ///////////////////////////////////////// Dahu ///////////////////////////////////////////////

    image2d<uint8_t> roi_alt_T = roi_alt(T, F, p_topleft, p_center);
    io::imsave(roi_alt_T, "roi_alt.pgm");

    image2d<uint8_t> roi_alt_new_T = roi_alt_new(T, F, p_topleft, p_center);
    io::imsave(roi_alt_new_T, "roi_alt_new.pgm");

    image2d<uint8_t> roi_T = roi(T, F, p_topleft, p_center);
    io::imsave(roi_T, "roi.pgm");



//    image2d<uint8_t> roi_fusion(roi_T.domain());

//    mln_foreach (auto p , roi_T.domain())
//    {
//        if (roi_T(p) == 255 and pathimage_dahu(p) == 255)
//        {
//            roi_fusion(p) =  255 ;
//        }
//        else
//        {
//            roi_fusion(p) =  0 ;
//        }
//    }

//    io::imsave(roi_fusion, "roi_fusion.pgm");








    //std::string fileout = outputDirectory + "/" + "fusion_" + fileName.substr(fileName.find("/") + 1);;
//    io::imsave(I, fileout.c_str());


////    fileName = outputDirectory + "/" + std::string(_dirent->d_name);
////    io::imsave(distance_map, fileName.c_str());

////    fileName = outputDirectory_scalar + "/" + std::string(_dirent->d_name);
////    io::imsave(distance_map_scalar, fileName.c_str());

        }
    }

    closedir(directory);

}

