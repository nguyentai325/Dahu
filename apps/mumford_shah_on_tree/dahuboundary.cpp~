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
    
    using namespace mln;
    typedef rgb8 V;
  
    using namespace mln::morpho::tos::internal;
    using namespace mln::morpho::tos;
    typedef irange<V> R;
  
    std::string inputDirectory = "/media/minh/DATA/Study/dataset/test_resized";
    std::string outputDirectory = "/home/minh/Documents/python/Evaluation/outpython/out_test_boundary";
    std::cout << inputDirectory << std::endl;
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

		if ((std::string(_dirent->d_name) != ".") and (std::string(_dirent->d_name) != ".."))
		{
			std::string fileName = inputDirectory + "/" +std::string(_dirent->d_name);
			std::cout << fileName << std::endl;
			//image2d<V> ima;
			//io::imread(fileName.c_str(), ima);	
			//cv::Mat ima = cv::imread(fileName.c_str());
			
			
			
			  //const char* input_path = argv[1];
			  int a0 = std::atoi(argv[2]);
			  int a1 = std::atoi(argv[3]);
			  int lambda = std::atoi(argv[4]);
			  const char* output_path = argv[5];


			  tbb::task_scheduler_init init;


			  // 1. Compute the individual ToS
			  

			  
			  image2d<V> ima;
			  io::imread(fileName.c_str(), ima);

			  image2d<V> f = addborder(ima, lexicographicalorder_less<value_t>());
			  //image2d<V> f = ima;
			  //io::imsave(f, "addborder.png");
			  
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
				//std::cout << "size T   "  << olddom1 << std::endl;
				
				
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
				//std::cout << "compute the Gos"  << std::endl;

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
				//std::cout << "size 2   "  << olddom2 << std::endl; 
			   
				
				std::ofstream file("abc.dot",
						 std::ios_base::out | std::ios_base::trunc);
				write_graphviz(file, T);    
			  }


				unsigned nodenumbers = T.nodes().size();
				
				
				const unsigned N1 = T.nodes().size();
				//std::cout << "N1  "  << N1 << std::endl; 

				auto& nodes1 = T._get_data()->m_nodes;
				auto& S     = T._get_data()->m_S;
				auto& pmap  = T._get_data()->m_pmap;	

				box2d D = pmap.domain();
				image2d<int> distancesum(D );
				image2d<uint8_t> distancefinal(D );
				
				extension::fill(distancesum, 0);
				

				
				//for (unsigned k = 0; k< 4 ;++k)
				//{
					//point2d startpoint= point2d(int(F.nrows()/2), int(F.ncols()/2));
					//if ((k==0))
					//	std::vector <point2d >   seedpoints(F.ncols());
					//if ((k==2))
					//	std::vector <point2d >   seedpoints(F.ncols());			
					//if ((k == 1))
					//	std::vector <point2d >   seedpoints(F.nrows());
					//if ((k==3))
					//	std::vector <point2d >   seedpoints(F.nrows());			
					
					std::queue<unsigned>  Q;
					int enqueue[nodenumbers] = {0};
					enqueue[0]  = 1;	
					int enqueue1[nodenumbers] = {0};
					enqueue1[0]  = 1;	    
					


					std::vector <point2d >   seedpoints(F.ncols()+ F.nrows());
					
					for (unsigned i = 0; i < F.ncols() ; ++i)
					{	
						
						seedpoints[i] = point2d(0, i);
						if ( enqueue[pmap(seedpoints[i])] == 0)
						{
							Q.push(pmap(seedpoints[i]));
							enqueue[pmap(seedpoints[i])] = 1;
							enqueue1[pmap(seedpoints[i])] = 1;
							
						}
					}
	 
					for (unsigned i = 0; i < F.nrows() ; ++i)
					{
						seedpoints[i] = point2d(i, 0);
						if ( enqueue[pmap(seedpoints[i])] == 0)
						{
							Q.push(pmap(seedpoints[i]));
							enqueue[pmap(seedpoints[i])] = 1;
							enqueue1[pmap(seedpoints[i])] = 1;
							
						}				
					}
					
					for (unsigned i = 0; i < F.ncols() ; ++i)
					{
						seedpoints[i] = point2d(F.nrows()-1, i);
						if ( enqueue[pmap(seedpoints[i])] == 0)
						{
							Q.push(pmap(seedpoints[i]));
							enqueue[pmap(seedpoints[i])] = 1;
							enqueue1[pmap(seedpoints[i])] = 1;
							
						}			
					}
					
					for (unsigned i = 0; i < F.nrows() ; ++i)
					{
						seedpoints[i] = point2d(i, F.ncols()-1);
						if ( enqueue[pmap(seedpoints[i])] == 0)
						{
							Q.push(pmap(seedpoints[i]));
							enqueue[pmap(seedpoints[i])] = 1;
							enqueue1[pmap(seedpoints[i])] = 1;
							
						}
					}																		

					
					//while(not Q1.empty())
					//{
					//	unsigned nodeid = Q1.front();
					//	Q1.pop();
					//	std::cout << "node id  "  << nodeid << std::endl;		
					//}
					
					

					//std::cout << "node numbers   "  << T.nodes().size() << std::endl;
					
					//std::cout <<  "start " <<pmap(startpoint)    << std::endl; 
					
					image2d<V> distancemap(D );
					image2d<V> min_pixel(D );
					image2d<V> max_pixel(D );
					
					
					image2d<uint8_t> distancemap1(D );
					
					//unsigned start = pmap(startpoint);
					
						


					

					
					std::vector<V>  min_node(nodenumbers);
					std::vector<V>  max_node(nodenumbers);
					std::vector<V>  dist(nodenumbers);	
					
					

					std::vector< std::vector<unsigned> > children(nodenumbers);
					
					//int enqueue[nodenumbers] = {0};
					//enqueue[0]  = 1;


					//std::cout << "initiation   "  << std::endl;

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
						//std::cout << "median value"  << medianvalue << "  " << medianvalue1  << std::endl;


						for (int i = 0; i< 3; ++i)
						{
							min_node[x.id()][i] = medianvalue1[i];
							max_node[x.id()][i] = medianvalue1[i];
							dist[x.id()][i] = max_node[x.id()][i] - min_node[x.id()][i];			
						}		
						
						unsigned q = x.get_parent_id();
						children[q].push_back(x.id());
					}


					//std::queue<unsigned>  Q;
					
					//Q.push(start);
					//enqueue[start] = 1;
					
					//auto x1 = T.get_node(start);
					//
					//std::cout << "x1.id  "  <<x1.id()  << std::endl; 
				  
				  
					//unsigned par_p = (T.get_node(start).get_parent_id());    
					//std::cout << "par_p   "  << par_p  << std::endl;
					
					//std::cout <<  "propagation "  << std::endl;

					while(not Q.empty())
					{
						unsigned nodeid = Q.front();
						Q.pop();

						for (int i = 0; i< 3; ++i)
						{
							dist[nodeid][i] = max_node[nodeid][i] - min_node[nodeid][i];			
						}	
						
						unsigned par_p = (T.get_node(nodeid).get_parent_id());

						if ((enqueue[par_p] == 0) and (enqueue1[par_p] == 1))
						{
							for (int i = 0; i< 3; ++i)
							{
								if (max_node[par_p][i] - min_node[par_p][i] > max_node[nodeid][i] - min_node[nodeid][i])
								{
									max_node[par_p][i] = max_node[nodeid][i];
									min_node[par_p][i] = min_node[nodeid][i];
									
								}
							}

						}

								
						if ((enqueue[par_p] == 0) and (enqueue1[par_p] == 0))
						{
							for (int i = 0; i< 3; ++i)
							{
								if (min_node[par_p][i] > min_node[nodeid][i])
									min_node[par_p][i] = min_node[nodeid][i];
								if (max_node[par_p][i] < max_node[nodeid][i])
									max_node[par_p][i] = max_node[nodeid][i];
							}
							Q.push(par_p);
							enqueue1[par_p] =1;
						}


						const unsigned nchildren = children[nodeid].size();
						
						for (unsigned j= 0; j< nchildren ; ++j)
						{
							unsigned child = children[nodeid][j];


							if ((enqueue[child] == 0) and (enqueue1[child] == 1) )
							{
								for (int i = 0; i< 3; ++i)
								{
									if (max_node[child][i] - min_node[child][i] > max_node[nodeid][i] - min_node[nodeid][i])
									{
										max_node[child][i] = max_node[nodeid][i];
										min_node[child][i] = min_node[nodeid][i];
										
									}
								}
						
							}

										
							if ((enqueue[child] == 0) and (enqueue1[child] == 0) )
							{
								for (int i = 0; i< 3; ++i)
								{
									if (min_node[child][i] > min_node[nodeid][i])
										min_node[child][i] = min_node[nodeid][i];
									if (max_node[child][i] < max_node[nodeid][i])
										max_node[child][i] = max_node[nodeid][i];
								}
								Q.push(child);
								enqueue1[child] =1;				
							}
						
						}
						
					}
					
					
					//  back propagation
					
					
					mln_foreach(auto x, T.nodes())
					{
						mln_foreach (auto p, x.proper_pset())
						{
							for (int i = 0; i< 3; ++i)
							{
								min_pixel(F.point_at_index(p))[i] = min_node[x.id()][i];
								max_pixel(F.point_at_index(p))[i] = max_node[x.id()][i];		
								distancemap(F.point_at_index(p))[i] = 	max_pixel(F.point_at_index(p))[i] - 	min_pixel(F.point_at_index(p))[i] ;				
							}
							distancemap1(F.point_at_index(p)) = distancemap(F.point_at_index(p))[0]/3 + distancemap(F.point_at_index(p))[1]/3+ distancemap(F.point_at_index(p))[2]/3;
							//distancesum(F.point_at_index(p)) = distancesum(F.point_at_index(p)) + distancemap1(F.point_at_index(p));
							//distancefinal(F.point_at_index(p)) = distancesum(F.point_at_index(p)) /4 ;
							
						}

					}	
					

					//io::imsave(distancemap, "dmap0.png");

					
		
			
			
			
			
			
			
			fileName = outputDirectory + "/" + std::string(_dirent->d_name);			
			
			io::imsave(distancemap1, fileName.c_str()); 		
			//cv::imwrite(fileName.c_str(), ima);
		}	
		
    }
    
    closedir(directory);    

	
}
