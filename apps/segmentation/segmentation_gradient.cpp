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
        if (p[0] == 0 or p[0]== height -1 or p[1] == 0 or p[1] == width)
            return true;
        else
            return false;

    }

    bool is_border1 (point2d p, unsigned height , unsigned width)
    {
        if (p[0] != 0  and p[0]  != 1  and p[0] != height -1 and p[0] != height -2 and p[1] != 0  and p[1]  != 1 and p[1]  != width - 1 and p[1] != width -2)
        {
            if (p[0] == 2 or p[0] == height - 3 or p[1] == width -3 or p[1] == 2  )
                return true;
            else
                return false;
        }
    }

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

    image2d<uint8> interpolate(const image2d<uint8>& ima)
    {
      image2d<uint8> out(2*ima.nrows()-1, 2*ima.ncols()-1);
      typedef point2d P;
      mln_foreach(point2d p, ima.domain())
        {
      uint8 a = ima.at(p),
        b = ima.at(p + P{0,1}),
        c = ima.at(p + P{1,0}),
        d = ima.at(p + P{1,1});

      point2d q = 2 * p;
      out.at(q) = ima.at(p);
      out.at(q + P{0,1}) = (a + b) / 2;
      out.at(q + P{1,0}) = (a + c) / 2;
      out.at(q + P{1,1}) = (a + b + c + d) / 4;
        }

      return out;
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






  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;

  image2d<V> ima;
  io::imread(input_path, ima);

  image2d<uint8>  ima_gray = rgb2gray(ima);






  image2d<uint8> f = addborder_gray(ima_gray);

  image2d<uint8> F = interpolate(f);


  tree_t T;
  T = morpho::cToS(f, c4);
  if (a0 > 0)
  {
    grain_filter_inplace(T, a0);
    T.shrink_to_fit();
  }
  T._reorder_pset();



  auto& U  = T._get_data()->m_Uv;

  // Initiation

  unsigned nodenumbers = T.nodes().size();
  std::cout << "node numbers "  << nodenumbers << std::endl;

  auto& nodes1 = T._get_data()->m_nodes;
  auto& S     = T._get_data()->m_S;
  auto& pmap  = T._get_data()->m_pmap;

  box2d D = pmap.domain();
  std::cout <<"Domain  " <<D << std::endl;
  unsigned height = D.pmax[0];
  unsigned width = D.pmax[1];



//  std::vector<int>  depth(nodenumbers);
//  std::fill (depth.begin(),depth.end(),0);


//  mln_foreach(auto x, T.nodes())
//  {
//      if (x.id()>1)
//      {
//        depth[x.id()] = depth[x.get_parent_id()] + 1;
//        std::cout << x.id()  << " has parent "  << x.get_parent_id()   << "  depth  " << depth[x.id()]  << std::endl;
//      }
//  }



  std::vector<int>  area(nodenumbers);
  std::fill (area.begin(),area.end(),0);

  std::vector<int>  average_gradient(nodenumbers);
  std::fill (average_gradient.begin(),average_gradient.end(),0);

//  accu::accumulators::accu_if<accu::accumulators::count<>,
//                              K1::is_face_2_t, point2d> acc;
//  auto area = morpho::paccumulate(T, T._get_data()->m_pmap, acc);

  std::vector<std::vector<point2d> >  gradient(nodenumbers);



  //std::cout << "index   "  << F.point_at_index(F.index_of_point(point2d(1,1)))  << std::endl;



  mln_foreach(auto x, T.nodes())
  {
      image2d <uint8>  node_ima(D);
      extension::fill(node_ima, 0);
//      area[x.id()] = x.proper_pset().size();
      int sum_gradient = 0;

      mln_foreach (auto p, x.proper_pset())
      {
          node_ima(F.point_at_index(p)) = 255;
          area[x.id()] = area[x.id()] + 1;
          const int dx4[4] = {-1,  0,  1,  0};
          const int dy4[4] = { 0, -1,  0,  1};
          point2d p1 = F.point_at_index(p);

          bool check = false;
          int grad = 0;
          for( int n = 0; n < 4; n++ )
          {
                int pt_y = p1[0] + dy4[n];
                int pt_x = p1[1] + dx4[n];
                if( (pt_x >= 0 && pt_x < width) && (pt_y >= 0 && pt_y < height) )
                {
                    point2d neibor = point2d(pt_y,pt_x);
                    int n1 = F.index_of_point(neibor);
                    if (T.get_node_id(n1) != x.id())
                    {
                        check = true;
                        int temp = std::abs(U(U.point_at_index(p)) - U(U.point_at_index(n1)));
                        if (grad < temp)
                            grad = temp;
                    }
                }
          }
          if (check == true)
          {
              gradient[x.id()].push_back(p);
              sum_gradient = sum_gradient + grad;
          }



      }
      if (area[x.id()] > 100 and area[x.id()]/gradient[x.id()].size() > 4)
          average_gradient[x.id()] = sum_gradient/ gradient[x.id()].size();
      else
          average_gradient[x.id()] = 0;

      std::cout << x.id()  << "   " << gradient[x.id()].size() <<"  " << area[x.id()] << "  "  << average_gradient[x.id()]<< std::endl;



      std::string fileName = std::to_string(x.id()) +  ".png";
      io::imsave(node_ima, fileName.c_str());
  }








}

