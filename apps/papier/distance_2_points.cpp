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



#include <sstream>
#include <stdlib.h>


#include <iostream>
#include <stdio.h>
#include "dirent.h"


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

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 or  p[0] == height -1 or p[1] == 0  or p[1] == width -1 )
            return true;
        else
            return false;
    }

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

    template <typename V>
    image2d<V> waterflow_mbd_distance_2_points(image2d<V> I, point2d p_topleft, point2d p_bottomright, V & waterflow_distance )
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
                if (i == p_topleft[0]/2 and j == p_topleft[1]/2 )
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
        waterflow_distance  = dmap(p_bottomright/2);

        return dmap;

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


  std::string inputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/test_resized";
  std::string outputDirectory = "/home/minh/Documents/code/code_edwin/build/apps/papier/distance2points";
  std::string outputDirectory_scalar = "/home/minh/Documents/code/code_edwin/build/apps/papier/output_dahu_color/scalar";

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

//  while((_dirent = readdir(directory)) != NULL)
//  {

//    if ((std::string(_dirent->d_name) != ".") and (std::string(_dirent->d_name) != "..") )
//    {


  std::string fileName = input_path;


  //std::string fileName = input_path;

  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<V> ima;
  io::imread(fileName.c_str(), ima);

  image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
  image2d<V> F = interpolate_k1(f);
  unsigned height_F = F.nrows();
  unsigned width_F = F.ncols();

  //std::string outputDirectory = "/lrde/home/movn/Documents/code/code_edwin2/program/apps/mumford_shah_on_tree/abc";



  //////////////////////////////////   Tree of shapes ///////////////////////////////////////////////////////////////

  image2d<V> ima_compo(F.domain() );



  // COmpute the color tree of shapes T
  tree_t T;
  {
    // 1. Compute the marginal ToS and filter them if necessary.
    tree_t t[NTREE];
    tbb::parallel_for(0, (int)NTREE, [&t,&f,a0](int i){
        t[i] = morpho::cToS(imtransform(f, [i](value_t x) { return x[i]; }), c4);
        if (a0 > 0) {
          grain_filter_inplace(t[i], a0);
          t[i].shrink_to_fit();
        }
      });


    auto& U0  = t[0]._get_data()->m_Uv;
    auto& U1  = t[1]._get_data()->m_Uv;
    auto& U2  = t[2]._get_data()->m_Uv;


    mln_foreach(point2d p1, U0.domain())
    {
        ima_compo(p1)[0] = U0(p1);
        ima_compo(p1)[1] = U1(p1);
        ima_compo(p1)[2] = U2(p1);
    }

    // 2. Compute the Gos.
    MyGraph g2;
    std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
    std::tie(g2, tlink) = compute_g2<NTREE>(t);

    // 3. Compute the depth image
    boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
    image2d<uint16> imdepth = imchvalue<uint16>(F);

    write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

    // debug
    // io::imsave(imdepth, "depth.tiff");

    // 4. Compute the saturated maxtree
    std::tie(T, std::ignore) = satmaxtree(imdepth);
    T = tree_keep_2F(T);
    T._reorder_pset();





    std::ofstream file("abc.dot",
             std::ios_base::out | std::ios_base::trunc);
    write_graphviz(file, T);

  }



    // 2 points

    point2d p_topleft = point2d(130,130);
    point2d p_bottomright = point2d(height_F-79,width_F-79);
    point2d p_center = point2d(272, 364);


    V dahu_distance;
    image2d<V> dahu_map =  dahu_distance_2_points( T,  ima_compo,  p_topleft,  p_bottomright, dahu_distance );
    io::imsave(dahu_map, "dmap1.png");

    std::cout << p_bottomright << std::endl;
    std::cout << int(dahu_map(p_bottomright)[0])  << "  " << int(dahu_map(p_bottomright)[1]) << "  "  <<int(dahu_map(p_bottomright)[2])  << std::endl;
    std::cout << int(dahu_map(p_center)[0])  << "  " << int(dahu_map(p_center)[1]) << "  "  <<int(dahu_map(p_center)[2])  << std::endl;


    ///////////////////////////////////////  Waterflow MBD  ////////////////////////////////////////////////////////////


    image2d<V> I = addborder_color(ima);

    V waterflow_distance;
    image2d<V> waterflow_map = waterflow_mbd_distance_2_points( I,  p_topleft,  p_bottomright, waterflow_distance );
    io::imsave(waterflow_map, "dmap2.png");

    std::cout << int(waterflow_map(p_bottomright/2)[0])  << "  " << int(waterflow_map(p_bottomright/2)[1]) << "  "  <<int(waterflow_map(p_bottomright/2)[2])  << std::endl;
    std::cout << int(waterflow_map(p_center/2)[0])  << "  " << int(waterflow_map(p_center/2)[1]) << "  "  <<int(waterflow_map(p_center/2)[2])  << std::endl;


//    I(point2d(50,50)) = {255,0,0};
//    I(point2d(201,284)) = {255,0,0};

    for (int i = 0; i < I.nrows(); i++)
    {
        for(int j = 0; j < I.ncols() ; j++)
        {
            point2d p = point2d(i,j);
            if ((i > p_topleft[0]/2 - 4 and i < p_topleft[0]/2 + 4 and j > p_topleft[1]/2 - 4 and j < p_topleft[1]/2 + 4)  or (i > p_bottomright[0]/2-4 and i < p_bottomright[0]/2+4 and j > p_bottomright[1]/2-4 and j < p_bottomright[1]/2+4) )
            {
                I(point2d(i,j)) = {255,0,0};
            }

            if ((i > p_center[0]/2-4 and i < p_center[0]/2+4 and j > p_center[1]/2-4 and j < p_center[1]/2+4)  )
            {
                I(point2d(i,j)) = {0,0,255};
            }

        }
    }

    io::imsave(I, output_path);



    ///////////////////////////////////////  Fast MBD  ////////////////////////////////////////////////////////////









    ///////////////////////////////////////  Fast Dahu  ////////////////////////////////////////////////////////////


    std::string fileout;
    fileout = outputDirectory + "/" + "dahu_" + fileName.substr(fileName.find("/") + 1);
    std::cout << fileout << std::endl;
    io::imsave(dahu_map, fileout.c_str());


    fileout = outputDirectory + "/" + "water_" + fileName.substr(fileName.find("/") + 1);;
    io::imsave(waterflow_map, fileout.c_str());


//    fileName = outputDirectory + "/" + std::string(_dirent->d_name);
//    io::imsave(distance_map, fileName.c_str());

//    fileName = outputDirectory_scalar + "/" + std::string(_dirent->d_name);
//    io::imsave(distance_map_scalar, fileName.c_str());

//        }
//    }

//    closedir(directory);

}

