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

  /////////////////////////////////////////
  for(int i = 0; i< 5; ++i)
    for (int j = 0; j<5 ; ++j)
       std::cout << "pixel "  << int(ima(point2d(i,j))[1])  << std::endl;
  
  
  
  ////////////////////////////////////////////
  

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
	
    io::imsave(ima_compo, "ima_compo.png");
    io::imsave(U0, "U0.png");
    io::imsave(U1, "U1.png");
    io::imsave(U2, "U2.png");
    
    
    
    

	
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
	
	
    t[0]._reorder_pset();
    t[1]._reorder_pset();
    t[2]._reorder_pset();


	
	for (int k = 0; k < 3; ++k)
	{
	
		const unsigned N1 = t[k].nodes().size();
		std::cout << "N1  "  << N1 << std::endl; 

		auto& nodes1 = t[k]._get_data()->m_nodes;
		auto& S     = t[k]._get_data()->m_S;
		auto& pmap  = t[k]._get_data()->m_pmap;	

		

		
		//box2d b = pmap.domain();
		int b = nodes1.size();
		//std::cout << "b   "  << b << std::endl;

		

		
		
		

		
		
		//  compute dahu distance
		
		// initiation
		
		box2d D = pmap.domain();
		unsigned N = S.size();
		image2d<uint8_t> distancemap(D );
		image2d<uint8_t> min_pixel(D );
		image2d<uint8_t> max_pixel(D );
		
		
		image2d<uint8_t> distancemap1(D );
		
		
		std::cout << " D "  << D << std::endl;
		std::cout << " N " << N << std::endl;
		
		// take started node
		
		unsigned start,pa;
		mln_foreach(auto x, t[k].nodes())
		{
			start = x.id();
		
			break;
		}    
		std::cout << "start   "  << start   << "pa "<< pa << std::endl;
		
		
		// number of nodes in the tree
		unsigned nodenumbers = 0;
			

		std::cout << "node numbers"  << t[k].nodes().size() << std::endl;
		nodenumbers = t[k].nodes().size();
		
		int see[nodenumbers] = {0};
		see[start] = 1;
		
		std::vector<uint8_t>  min_node(nodenumbers);
		std::vector<uint8_t>  max_node(nodenumbers);
		std::vector<uint8_t>  dist(nodenumbers);
		
		

		// set the initial value to each node 
		
		mln_foreach(auto x, t[k].nodes())
		{
			
			std::vector<uint8_t> listpixel1;
			
			mln_foreach (auto p, x.proper_pset())
			{			
				listpixel1.push_back(ima_compo(U2.point_at_index(p))[k]);
				
				//std::cout << "check " << x.id() << "  "<<F.point_at_index(p)  << "  "  << ima_compo(F.point_at_index(p)) << std::endl;
			}
		
			
			std::partial_sort(listpixel1.begin(), listpixel1.begin() + listpixel1.size()/2+1, listpixel1.end());

			uint8_t medianvalue1;

			medianvalue1 = listpixel1[listpixel1.size()/2];
			//std::cout << "median value"  << medianvalue << "  " << medianvalue1  << std::endl;


			min_node[x.id()] = medianvalue1;
			max_node[x.id()] = medianvalue1;
			dist[x.id()] = max_node[x.id()] - min_node[x.id()];			
				
		
		}	
		
		
		
		

		// update value of each node by comparing with its parent
		 
		std::cout << "propagate "  << std::endl;


		

		mln_foreach(auto x, t[k].nodes())
		{
			if (see[x.id()] == 0)
			{
				if (min_node[x.id()] > min_node[x.get_parent_id()])
					min_node[x.id()] = min_node[x.get_parent_id()];
				if (max_node[x.id()] < max_node[x.get_parent_id()])
					max_node[x.id()] = max_node[x.get_parent_id()];
			}
		}	
		
		
		// propagate dahu distance to every pixel in the image
		
		mln_foreach(auto x, t[k].nodes())
		{
			mln_foreach (auto p, x.proper_pset())
			{
				min_pixel(F.point_at_index(p)) = min_node[x.id()];
				max_pixel(F.point_at_index(p)) = max_node[x.id()];		
				distancemap(F.point_at_index(p))= 	max_pixel(F.point_at_index(p)) - 	min_pixel(F.point_at_index(p)) ;				

				
			}

		}			
		if (k== 0)
			io::imsave(distancemap, "dmap_test0.png");
		if (k== 1)
			io::imsave(distancemap, "dmap_test1.png");
		if (k== 2)
			io::imsave(distancemap, "dmap_test2.png");			
	}
    
    
    
    


	///////////////////////////////////////////////////////////////////

    // 2. Compute the Gos.
    MyGraph g2;
    std::array<property_map<tree_t, typename MyGraph::vertex_descriptor>, NTREE> tlink;
    std::tie(g2, tlink) = compute_g2<NTREE>(t);

    // 3. Compute the depth image
    boost::vector_property_map<unsigned> gdepth = compute_graph_depth(g2);
    image2d<uint16> imdepth = imchvalue<uint16>(F);

	unsigned sz10 = imdepth.domain().size();
    std::cout <<  "sz10   "   << sz10   << std::endl;	

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
	write_graphviz(file, t[1]);    
  }



	const unsigned N1 = T.nodes().size();
	std::cout << "N1  "  << N1 << std::endl; 

	auto& nodes1 = T._get_data()->m_nodes;
	auto& S     = T._get_data()->m_S;
	auto& pmap  = T._get_data()->m_pmap;	

	

	
	//box2d b = pmap.domain();
	int b = nodes1.size();
	//std::cout << "b   "  << b << std::endl;

	

	
	
	

	
	
	//  compute minimum barier
	
	
	
	box2d D = pmap.domain();
	unsigned N = S.size();
	image2d<V> distancemap(D );
	image2d<V> min_pixel(D );
	image2d<V> max_pixel(D );
	
	
	image2d<uint8_t> distancemap1(D );
	//image2d<V> distancemap(D );
	
	
	std::cout << " D "  << D << std::endl;
	std::cout << " N " << N << std::endl;
	
    

    
    
    
    unsigned nodenumbers = 0;
        

    std::cout << "node numbers"  << T.nodes().size() << std::endl;
    nodenumbers = T.nodes().size();
    


	unsigned root = T.get_root_id();
	
	std::cout << "root "  << root <<  std::endl; 
	
	
	point2d startpoint  =  point2d(2032, 2668);
	unsigned start = pmap(startpoint);
	
	std::vector <bool>  under(nodenumbers, 0);
	std::fill(under.begin(), under.end(), false);	
	
	under[start] = true;	

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
	
	io::imsave(under_ima, "cutimage.png");

    

    




}
