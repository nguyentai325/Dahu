#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/morpho/component_tree/io.hpp>

#include <mln/morpho/tos/ctos.hpp>
#include <mln/morpho/component_tree/compute_depth.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/graphviz.hpp>
#include <mln/morpho/extinction.hpp>

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
#include "mumford_shah_on_tree.hpp"
#include <iostream>
#include <fstream>

# include <mln/core/extension/fill.hpp>

# include <mln/core/image/image2d.hpp>
# include <mln/core/neighb2d.hpp>
# include <vector>
# include <mln/core/extension/fill.hpp>

# include <mln/morpho/pqueue_fast.hpp>

# include <mln/morpho/tos/pset.hpp>
# include <mln/morpho/tos/pset_priority.hpp>
# include <queue>

# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/fill.hpp>

# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/pqueue_fast.hpp>

# include <mln/core/wrt_offset.hpp>

# include <vector>
# include <stack>
# include <queue>




# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/algorithm/fill.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>
# include <mln/morpho/tos/pset.hpp>
# include <mln/morpho/tos/pset_priority.hpp>
# include <mln/morpho/maxtree/maxtree.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/core/image/morphers/casted_image.hpp>


# include <mln/core/image/image2d.hpp>
# include <mln/core/neighb2d.hpp>
# include <vector>
# include <mln/morpho/structural/dilate.hpp>




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

int increasing(int a, int b)
{
	if (a < b)
		return -1;
	return a > b;
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
 // int a1 = std::atoi(argv[3]);
  //int lambda = std::atoi(argv[4]);
  //const char* output_path = argv[5];


  tbb::task_scheduler_init init;


  // 1. Compute the individual ToS
  using namespace mln;
  typedef rgb8 V;
  
  using namespace mln::morpho::tos::internal;
  using namespace mln::morpho::tos;
  typedef irange<V> R;
  
  image2d<V> ima;
  io::imread(input_path, ima);

  image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
  //image2d<V> f = ima;
  io::imsave(f, "addborder.png");
  
  image2d<V> F = interpolate_k2(ima);


  

  image2d<V> ima_compo(F.domain() );
  
  

  
  
  

  // COmpute the color tree of shapes T
  tree_t T;
  {
    // 1. Compute the marginal ToS and filter them if necessary.
    tree_t t[NTREE];
    tbb::parallel_for(0, (int)NTREE, [&t,&ima,a0](int i){
        t[i] = morpho::cToS(imtransform(ima, [i](value_t x) { return x[i]; }), c4);
        if (a0 > 0) {
          grain_filter_inplace(t[i], a0);
          t[i].shrink_to_fit();
          

        }
      });
    
    auto newdata1 = t[1]._get_data();
    box2d olddom1 = newdata1->m_pmap.domain();
    std::cout << "size T   "  << olddom1 << std::endl;
    
    
    auto& U0  = t[0]._get_data()->m_Uv;
    auto& U1  = t[1]._get_data()->m_Uv;
    auto& U2  = t[2]._get_data()->m_Uv;
    
       
	mln_foreach(point2d p1, U0.domain())
	{
		ima_compo(p1)[0] = U0(p1);
		ima_compo(p1)[1] = U1(p1);
		ima_compo(p1)[2] = U2(p1);
	}
    

	
	///////////////////////////////////////////////////////////////
	
	// Compute the dahu distance in each color of original image 
	
	// 1 node <=>  set of pixels
	// choose 1 pixel to represent for all the pixel in the node
	// take its value of sorting step
	// computing dahu distance 
	// create a loop to touch every node in the tree
	// update max and min value in each node
	// dahu distance = max - min
	// propagate value of dahu distance to every pixel in the node 
	
	

    


	///////////////////////////////////////////////////////////////////

    // 2. Compute the Gos.
    MyGraph g2;
    std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
    std::tie(g2, tlink) = compute_g2<NTREE>(t);

    // 3. Compute the depth image
    boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
    image2d<uint16> imdepth = imchvalue<uint16>(F);

	unsigned sz10 = imdepth.domain().size();

    write_vmap_to_image(t, &tlink[0], gdepth, imdepth);

    // debug

    //io::imsave(imdepth, "depth.tiff");

    // 4. Compute the saturated maxtree
    std::tie(T, std::ignore) = satmaxtree(imdepth);
	
    

    
    
    T = tree_keep_2F(T);
    T._reorder_pset();
 
    auto newdata2 = T._get_data();
    auto olddom2 = newdata2->m_S.size();
    std::cout << "size 2   "  << olddom2 << std::endl; 
   
    
    std::ofstream file("abc.dot",
			 std::ios_base::out | std::ios_base::trunc);
	write_graphviz(file, T);    
  }

   
	const unsigned N1 = T.nodes().size();
	std::cout << "N1  "  << N1 << std::endl; 

	auto& nodes1 = T._get_data()->m_nodes;
	auto& S     = T._get_data()->m_S;
	auto& pmap  = T._get_data()->m_pmap;
	
	
	box2d D = pmap.domain();
	unsigned nodenumbers = T.nodes().size();   
   
    std::cout << "node number  "  << nodenumbers <<  std::endl;
    // Compute depth image
  
	image2d <unsigned>  depth(D);	
	std::vector<unsigned> depthvec(nodenumbers, 0);
	std::fill(depthvec.begin(), depthvec.end(), 0);
	
	extension::fill(depth, 0);
	
	unsigned root = T.get_root_id();
	
	std::cout << "root "  << root <<  std::endl; 
	

	
   	mln_foreach(auto x, T.nodes())
	{
		if (x.id()  == root)
		{
			mln_foreach (auto p, x.proper_pset())
			{
				depth(F.point_at_index(p)) = 0;
			}			
		}
		
		if (x.id() != root)
		{
			depthvec[x.id()]  = depthvec[x.get_parent_id()]  + 1;

			mln_foreach (auto p, x.proper_pset())
			{
				depth(F.point_at_index(p)) = depthvec[x.id()];
			}
		}
	}
  
  
    //mln_foreach(auto p , D)
    //{
		//std::cout << "p  " << p   << "   "  << depth(p) <<  "    "  << pmap(p)  << std::endl;  

	//}



	// compute_tree_path
	//
	image2d <uint8_t>  out_path(D);
	extension::fill(out_path, 0);		
	
	for (unsigned k= 0; k < F.ncols()-1 ; ++k)
	{
		std::cout << "k  "  << k  << std::endl;
		point2d startpoint  =  point2d(0, k);
		point2d endpoint    =  point2d(F.nrows()-1, k);
		//std::cout << "startpoint   "  <<  startpoint   << "    "  << endpoint  << std::endl; 
		unsigned lca;
		
		
		std::vector <unsigned >  path_;

		unsigned start = pmap(startpoint);
		unsigned end = pmap(endpoint);
		
		//std::cout << "start  "  << start << "end  "  << end  << std::endl;
		
		unsigned
			p_cur = depthvec[start] >= depthvec[end]   ? start :   end,
			p_obj = depthvec[start] >= depthvec[end] ? end   : start;	
		//std::cout << "p_cur  "  << p_cur << "p_obj  "  << p_obj  << std::endl;
			
		unsigned
			d_cur = depthvec[p_cur],
			d_obj = depthvec[p_obj]; 
		//std::cout << "d_cur  "  << d_cur << "d_obj  "  << d_obj  << std::endl;
			 

		while (d_cur > d_obj)
		{
			d_cur -= 1;
			path_.push_back(p_cur);
			p_cur = T.get_node(p_cur).get_parent_id();
		}    	
			  
			
		while (p_cur != p_obj)
		{
			path_.push_back(p_cur);
			path_.push_back(p_obj);
			p_cur = T.get_node(p_cur).get_parent_id();
			p_obj = T.get_node(p_obj).get_parent_id();
		}

		path_.push_back(p_cur);
		lca = p_cur;	
		//std:: cout << "lca  "  << lca << std::endl;
		
		
		
		// show tree path 
		
		
		std::vector<bool> path_ima(nodenumbers, 0);
		image2d<bool> imagepath(D);
		
		std::fill(path_ima.begin(), path_ima.end(), false);	
		
		
		unsigned pathsize = path_.size();
		for (unsigned i = 0; i < pathsize ; ++i)
		{
			//std::cout << "path   "  << path_[i]  << std::endl;
			path_ima[path_[i]] = true;
		}	

		std::vector <bool>  under(nodenumbers, 0);
		std::fill(under.begin(), under.end(), false);	
		
		mln_foreach(auto x, T.nodes())
		{
			if ((path_ima[x.id()] == false) and (path_ima[x.get_parent_id()]== true))
			{
				under[x.id()] = true;
				//std::cout << "under  "  << x.id() << "   "  << under[x.id()] << std::endl;
			}
		}
		
		
		
		// back propagate  subtree (under)
		image2d <bool > under_ima(D);
		extension::fill(under_ima, false);
		
		
		mln_foreach(auto x, T.nodes())
		{
			if ((under[x.id()] == false) and (under[x.get_parent_id()] == true))
			{
				under[x.id()] = under[x.get_parent_id()];
				
				mln_foreach (auto p, x.proper_pset())
				{
					under_ima(F.point_at_index(p)) = under[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << under_ima(F.point_at_index(p)) << std::endl;
					
				}
			}
			
			else
			{
				mln_foreach (auto p, x.proper_pset())
				{
					under_ima(F.point_at_index(p)) = under[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << under_ima(F.point_at_index(p)) << std::endl;
					
				}			
			}
		}
		
		io::imsave(under_ima, "dmap_path.png");
		
		
		// back propagate subtree ( res)
		image2d <bool>  res_ima(D);
		extension::fill(res_ima, false);
		
		std::vector <bool>  res(nodenumbers, 0);
		std::fill(res.begin(), res.end(), false);	
		
		res[lca] = true;
		
		
		mln_foreach(auto x, T.nodes())
		{
			if ((res[x.id()] == false) and (res[x.get_parent_id()] == true))
			{
				res[x.id()] = res[x.get_parent_id()];
				
				mln_foreach (auto p, x.proper_pset())
				{
					res_ima(F.point_at_index(p)) = res[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << res_ima(F.point_at_index(p)) << std::endl;
					
				}
			}
			
			else
			{
				mln_foreach (auto p, x.proper_pset())
				{
					res_ima(F.point_at_index(p)) = res[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << res_ima(F.point_at_index(p)) << std::endl;
					
				}			
			}
		}	
		
		rect2d r = make_rectangle2d(3, 3);
		res_ima = morpho::structural::dilate(res_ima, r);
			 
		
		io::imsave(res_ima, "dmap_path1.png");

		image2d <bool>  result_ima(D);
		extension::fill(result_ima, false);
		
		mln_foreach (auto p , D)
		{
			if (res_ima(p) != under_ima(p))
			{
				result_ima(p) =  true ;
			}
		}

		io::imsave(result_ima, "dmap_path2.png");
		
		

		
		
		
		//   shortest path 
		
		image2d <int>  distancemap(D);
		extension::fill(distancemap, 0);
		//std::cout << "distancemap   " << distancemap(point2d(0,6))  << std::endl;
		
		


		image2d<bool>  dejavu(D);    
		extension::fill(dejavu, false); 
		
		std::queue<point2d>  Q;
		
		Q.push(startpoint);
		distancemap(startpoint)  = 1;
		dejavu(startpoint)  = true;
		
		//std::cout << "start point"  << startpoint[0]  << "   "  << startpoint[1]  << std::endl;
		
		// propagation
		
		while(not Q.empty())
		{
			point2d p1  = Q.front();
			Q.pop();
			if (p1 == endpoint)
			{
				break;
			}
			else
			{
				std::vector<point2d> neighbor(4, 0);
				neighbor[0] = point2d (p1[0]-1, p1[1]);
				neighbor[1] = point2d (p1[0]+1, p1[1]);
				neighbor[2] = point2d (p1[0], p1[1]-1);
				neighbor[3] = point2d (p1[0], p1[1]+1);
				mln_foreach(auto n1, neighbor)
				{
					if (D.has(n1) and dejavu(n1) == false   and result_ima(n1) == true)
					{
						//std::cout << "n1  "  << n1  << std::endl;
						distancemap(n1)  = distancemap(p1)  +1 ;
						dejavu(n1)  = true;
						Q.push(n1);
					}
				}

			}
			
		}
		
		//std::cout << "distancemap   " << distancemap(point2d(0,7))  << std::endl;
		//std::cout << "distancemap   " << distancemap(point2d(0,8))  << std::endl;
		//std::cout << "distancemap   " << distancemap(point2d(0,6))  << std::endl;

		// back trace
		//std::cout << "back trace"  << std::endl;
		
		std::vector<point2d>  minpath;
		point2d p2 = endpoint;
		out_path(p2) = out_path(p2) +1;
		while(p2 != startpoint)
		{
			std::vector <point2d> neighbor(4,0);
			neighbor[0] = point2d (p2[0]-1, p2[1]);
			neighbor[1] = point2d (p2[0]+1, p2[1]);
			neighbor[2] = point2d (p2[0], p2[1]-1);
			neighbor[3] = point2d (p2[0], p2[1]+1);	
			mln_foreach (auto n2, neighbor)
			{
				if (D.has(n2) and distancemap(p2) - distancemap(n2)== 1)
				{
					minpath.push_back(n2);
					//std::cout << "n2  " << n2  << std::endl; 
					p2 = n2;
					out_path(p2) = out_path(p2) + 1;
					break;
				}			
			}	
		}
	
	}
	
	
	for (unsigned k= 0; k < F.nrows()-1 ; ++k)
	{
		std::cout << "k  "  << k  << std::endl;
		point2d startpoint  =  point2d(k, 0);
		point2d endpoint    =  point2d(k, F.ncols()-1);
		//std::cout << "startpoint   "  <<  startpoint   << "    "  << endpoint  << std::endl; 
		unsigned lca;
		
		
		std::vector <unsigned >  path_;

		unsigned start = pmap(startpoint);
		unsigned end = pmap(endpoint);
		
		//std::cout << "start  "  << start << "end  "  << end  << std::endl;
		
		unsigned
			p_cur = depthvec[start] >= depthvec[end]   ? start :   end,
			p_obj = depthvec[start] >= depthvec[end] ? end   : start;	
		//std::cout << "p_cur  "  << p_cur << "p_obj  "  << p_obj  << std::endl;
			
		unsigned
			d_cur = depthvec[p_cur],
			d_obj = depthvec[p_obj]; 
		//std::cout << "d_cur  "  << d_cur << "d_obj  "  << d_obj  << std::endl;
			 

		while (d_cur > d_obj)
		{
			d_cur -= 1;
			path_.push_back(p_cur);
			p_cur = T.get_node(p_cur).get_parent_id();
		}    	
			  
			
		while (p_cur != p_obj)
		{
			path_.push_back(p_cur);
			path_.push_back(p_obj);
			p_cur = T.get_node(p_cur).get_parent_id();
			p_obj = T.get_node(p_obj).get_parent_id();
		}

		path_.push_back(p_cur);
		lca = p_cur;	
		//std:: cout << "lca  "  << lca << std::endl;
		
		
		
		// show tree path 
		
		
		std::vector<bool> path_ima(nodenumbers, 0);
		image2d<bool> imagepath(D);
		
		std::fill(path_ima.begin(), path_ima.end(), false);	
		
		
		unsigned pathsize = path_.size();
		for (unsigned i = 0; i < pathsize ; ++i)
		{
			//std::cout << "path   "  << path_[i]  << std::endl;
			path_ima[path_[i]] = true;
		}	

		std::vector <bool>  under(nodenumbers, 0);
		std::fill(under.begin(), under.end(), false);	
		
		mln_foreach(auto x, T.nodes())
		{
			if ((path_ima[x.id()] == false) and (path_ima[x.get_parent_id()]== true))
			{
				under[x.id()] = true;
				//std::cout << "under  "  << x.id() << "   "  << under[x.id()] << std::endl;
			}
		}
		
		
		
		// back propagate  subtree (under)
		image2d <bool > under_ima(D);
		extension::fill(under_ima, false);
		
		
		mln_foreach(auto x, T.nodes())
		{
			if ((under[x.id()] == false) and (under[x.get_parent_id()] == true))
			{
				under[x.id()] = under[x.get_parent_id()];
				
				mln_foreach (auto p, x.proper_pset())
				{
					under_ima(F.point_at_index(p)) = under[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << under_ima(F.point_at_index(p)) << std::endl;
					
				}
			}
			
			else
			{
				mln_foreach (auto p, x.proper_pset())
				{
					under_ima(F.point_at_index(p)) = under[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << under_ima(F.point_at_index(p)) << std::endl;
					
				}			
			}
		}
		
		io::imsave(under_ima, "dmap_path.png");
		
		
		// back propagate subtree ( res)
		image2d <bool>  res_ima(D);
		extension::fill(res_ima, false);
		
		std::vector <bool>  res(nodenumbers, 0);
		std::fill(res.begin(), res.end(), false);	
		
		res[lca] = true;
		
		
		mln_foreach(auto x, T.nodes())
		{
			if ((res[x.id()] == false) and (res[x.get_parent_id()] == true))
			{
				res[x.id()] = res[x.get_parent_id()];
				
				mln_foreach (auto p, x.proper_pset())
				{
					res_ima(F.point_at_index(p)) = res[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << res_ima(F.point_at_index(p)) << std::endl;
					
				}
			}
			
			else
			{
				mln_foreach (auto p, x.proper_pset())
				{
					res_ima(F.point_at_index(p)) = res[x.id()];
					//std::cout << F.point_at_index(p)  <<  "   "  << x.id() << "   "  << res_ima(F.point_at_index(p)) << std::endl;
					
				}			
			}
		}	
		
		rect2d r = make_rectangle2d(3, 3);
		res_ima = morpho::structural::dilate(res_ima, r);
			 
		
		io::imsave(res_ima, "dmap_path1.png");

		image2d <bool>  result_ima(D);
		extension::fill(result_ima, false);
		
		mln_foreach (auto p , D)
		{
			if (res_ima(p) != under_ima(p))
			{
				result_ima(p) =  true ;
			}
		}

		io::imsave(result_ima, "dmap_path2.png");
		
		

		
		
		
		//   shortest path 
		
		image2d <int>  distancemap(D);
		extension::fill(distancemap, 0);
		//std::cout << "distancemap   " << distancemap(point2d(0,6))  << std::endl;
		
		


		image2d<bool>  dejavu(D);    
		extension::fill(dejavu, false); 
		
		std::queue<point2d>  Q;
		
		Q.push(startpoint);
		distancemap(startpoint)  = 1;
		dejavu(startpoint)  = true;
		
		//std::cout << "start point"  << startpoint[0]  << "   "  << startpoint[1]  << std::endl;
		
		// propagation
		
		while(not Q.empty())
		{
			point2d p1  = Q.front();
			Q.pop();
			if (p1 == endpoint)
			{
				break;
			}
			else
			{
				std::vector<point2d> neighbor(4, 0);
				neighbor[0] = point2d (p1[0]-1, p1[1]);
				neighbor[1] = point2d (p1[0]+1, p1[1]);
				neighbor[2] = point2d (p1[0], p1[1]-1);
				neighbor[3] = point2d (p1[0], p1[1]+1);
				mln_foreach(auto n1, neighbor)
				{
					if (D.has(n1) and dejavu(n1) == false   and result_ima(n1) == true)
					{
						//std::cout << "n1  "  << n1  << std::endl;
						distancemap(n1)  = distancemap(p1)  +1 ;
						dejavu(n1)  = true;
						Q.push(n1);
					}
				}

			}
			
		}
		
		//std::cout << "distancemap   " << distancemap(point2d(0,7))  << std::endl;
		//std::cout << "distancemap   " << distancemap(point2d(0,8))  << std::endl;
		//std::cout << "distancemap   " << distancemap(point2d(0,6))  << std::endl;

		// back trace
		//std::cout << "back trace"  << std::endl;
		
		std::vector<point2d>  minpath;
		point2d p2 = endpoint;
		out_path(p2) = out_path(p2) +1;
		while(p2 != startpoint)
		{
			std::vector <point2d> neighbor(4,0);
			neighbor[0] = point2d (p2[0]-1, p2[1]);
			neighbor[1] = point2d (p2[0]+1, p2[1]);
			neighbor[2] = point2d (p2[0], p2[1]-1);
			neighbor[3] = point2d (p2[0], p2[1]+1);	
			mln_foreach (auto n2, neighbor)
			{
				if (D.has(n2) and distancemap(p2) - distancemap(n2)== 1)
				{
					minpath.push_back(n2);
					//std::cout << "n2  " << n2  << std::endl; 
					p2 = n2;
					out_path(p2) = out_path(p2) + 1;
					break;
				}			
			}	
		}
	
	}	
	io::imsave(out_path, "minpath.png");
	
	

}
