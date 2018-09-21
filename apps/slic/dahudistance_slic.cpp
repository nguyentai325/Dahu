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






namespace mln
{

    class Lab
    {
        public:

        Lab() {}
        Lab(float L, float a, float b) : L_(L), a_(a), b_(b) {}
        float  L() const { return L_; }
        float& L() { return L_; }
        float  a() const { return a_; }
        float& a() { return a_; }
        float  b() const { return b_; }
        float& b() { return b_; }
        private:
        float L_, a_, b_;
    };

    void rgb8_to_lab(const rgb8& c, float& L, float& a, float& b)
    {
        float
        R = float(c[0])   / 255.0,
        G = float(c[1]) / 255.0,
        B = float(c[2])  / 255.0;

        float X, Y,Z;
        float r1,g1,b1;
        float epsilon = 0.008856;	//actual CIE standard
        float kappa   = 903.3;		//actual CIE standard
        float Xr = 0.950456;	//reference white
        float Yr = 1.0;		//reference white
        float Zr = 1.088754;	//reference white
        double xr,yr,zr;
        double fx, fy, fz;

        if(R <= 0.04045)	r1 = R/12.92;
        else				r1 = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g1 = G/12.92;
        else				g1 = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b1 = B/12.92;
        else				b1 = pow((B+0.055)/1.055,2.4);


        X = r1*0.4124564 + g1*0.3575761 + b1*0.1804375;
        Y = r1*0.2126729 + g1*0.7151522 + b1*0.0721750;
        Z = r1*0.0193339 + g1*0.1191920 + b1*0.9503041;

        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;


        L = 116.0*fy-16.0;
        a = 500.0*(fx-fy);
        b = 200.0*(fy-fz);

    }

    void
    split_Lab(const image2d<Lab>& input,
            image2d<float>& L, image2d<float>& a, image2d<float>& b)
    {
        box2d D = input.domain();
        image2d<float> L_(D), a_(D), b_(D);

        mln_foreach(auto p, D)
        {
            L_(p) = input(p).L();
            a_(p) = input(p).a();
            b_(p) = input(p).b();
        }
        L = L_;
        a = a_;
        b = b_;
    }



    image2d<Lab>
    merge_Lab(const image2d<float>& L, const image2d<float>& a, const image2d<float>& b)
    {
        box2d D = L.domain();
        image2d<Lab> output(D);

        mln_foreach(auto p, D)
        {
            output(p).L() = L(p);
            output(p).a() = a(p);
            output(p).b() = b(p);
        }
        return output;
    }

    image2d<Lab> convert_and_shrink(const image2d<rgb8>& input)
    {
        image2d<Lab> output(input.nrows() ,
                            input.ncols() );
        box2d dom = output.domain();


        mln_foreach(auto p,dom)
        {
            float min = 1e10, max = -1e10, sum = 0;

            float L, a, b;

            rgb8_to_lab(input(p), L, a, b);
            if (L > max)
              max = L;
            if (a < min)
              min = a;
            sum += b;

            output(p).L() = max;
            output(p).a() = min;
            output(p).b() = sum;
        }

        return output;
    }


}




void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb] num  α₀ α₁ λ output[rgb]\n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}




int main(int argc, char** argv)
{
    if (argc < 6)
      usage(argv);

    double start_s=clock();

    const char* input_path = argv[1];
    int numSuperpixels = std::atoi(argv[2]);
    int a0 = std::atoi(argv[3]);
    int a2 = std::atoi(argv[4]);
    int lambda = std::atoi(argv[5]);
    const char* output_path = argv[6];
    std::string fileName = input_path;


    using namespace mln;
    typedef rgb8 V;

    image2d<V> Img;
    io::imread(fileName.c_str(), Img);

    //int numSuperpixels = 20;//default value
    double compactness = 20;//default value

    // convert rgb to Lab  *********************************


    image2d<Lab> input_Lab = convert_and_shrink(Img);
    image2d<float> L, a, b;
    split_Lab(input_Lab, L, a, b);

    box2d D = Img.domain();


    unsigned height = a.nrows();
    unsigned width = a.ncols();

    std::cout << height<< std::endl;
    std::cout << width << std::endl;

    int sz = height * width;
    int i,j;
    int index;
    float* lvec; float* avec; float* bvec;
    lvec    = (float*) malloc (sz*4);
    avec    = (float*) malloc (sz*4);
    bvec    = (float*) malloc (sz*4);

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width ; j++)
        {
            index = i*width + j;
            lvec[index] = L(point2d(i,j));
            avec[index] = a(point2d(i,j));
            bvec[index] = b(point2d(i,j));
        }
    }


    //---------------------------
    // Find seeds
    //---------------------------


    int STEP = sqrt((float)(sz)/(float)(numSuperpixels))+0.5;
    int n;
    int xstrips, ystrips;
    int xerr, yerr;
    float xerrperstrip,yerrperstrip;
    int xoff,yoff;
    int x,y;
    int xe,ye;
    int seedx,seedy;
    int i1;
    int* seedIndices;


    xstrips = (0.5+(float)(width)/(float)(STEP));
    ystrips = (0.5+(float)(height)/(float)(STEP));


    xerr = width  - STEP*xstrips;
    if(xerr < 0)
        {xstrips--;xerr = width - STEP*xstrips;}
    yerr = height - STEP*ystrips;
    if(yerr < 0)
        {ystrips--;yerr = height- STEP*ystrips;}

    xerrperstrip = (float)(xerr)/(float)(xstrips);
    yerrperstrip = (float)(yerr)/(float)(ystrips);

    xoff = STEP/2;
    yoff = STEP/2;

    seedIndices    = (int*) malloc (sz);

    n = 0;
    for( y = 0; y < ystrips; y++ )
    {
        ye = y*yerrperstrip;
        for( x = 0; x < xstrips; x++ )
        {
            xe = x*xerrperstrip;

            seedx = (x*STEP+xoff+xe);
            seedy = (y*STEP+yoff+ye);

            i1 = seedy*width + seedx;
            seedIndices[n] = i1;
            n++;
        }
    }
    int numseeds = n;



    float* kseedsx;float* kseedsy;
    float* kseedsl;float* kseedsa;float* kseedsb;

    kseedsx    = (float*) malloc (numseeds*4);
    kseedsy    = (float*) malloc (numseeds*4);
    kseedsl    = (float*) malloc (numseeds*4);
    kseedsa    = (float*) malloc (numseeds*4);
    kseedsb    = (float*) malloc (numseeds*4);

    int k;
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        kseedsl[k] = lvec[seedIndices[k]];
        kseedsa[k] = avec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];

    }


    //----------------------------------------------------------
    // Compute superpixels
    //-----------------------------------------------------------

    int x1, y1, x2, y2;
    float l1, a1, b1;
    float dist;
    float distxy;
    int itr;
    int n2;
    int x3,y3;
    int i2;
    int t;
    int index2;
    int ind;
    int r,c;
    int k1;
    int numk = numseeds;
    int offset = STEP;
    int* klabels;
    klabels = (int*) malloc (sz*4);//original k-means labels
    float* clustersize = (float*) malloc (numseeds*4);
    float* inv         = (float*) malloc (numseeds*4);
    float* sigmal      = (float*) malloc (numseeds*4);
    float* sigmaa      = (float*) malloc (numseeds*4);
    float* sigmab      = (float*) malloc (numseeds*4);
    float* sigmax      = (float*) malloc (numseeds*4);
    float* sigmay      = (float*) malloc (numseeds*4);
    float* distvec     = (float*) malloc (sz*4);




    float invwt = 1.0/((STEP/compactness)*(STEP/compactness));
    std::cout << " invwt "<< invwt << std::endl;



    for( itr = 0; itr < 10; itr++ )
    {
        for(t = 0; t < sz; t++)
        {distvec[t] = std::numeric_limits<float>::max();}

        for( n2 = 0; n2 < numk; n2++ )
        {
            x1 = kseedsx[n2]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n2]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n2]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n2]+offset; if(y2 > height) y2 = height;

            for( y3 = y1; y3 < y2; y3++ )
            {
                for( x3 = x1; x3 < x2; x3++ )
                {
                    index2 = y3*width + x3;

                    l1 = lvec[index2];
                    a1 = avec[index2];
                    b1 = bvec[index2];

                    dist =(l1 - kseedsl[n2])*(l1 - kseedsl[n2]) + (a1 - kseedsa[n2])*(a1 - kseedsa[n2]) +  (b1 - kseedsb[n2])*(b1 - kseedsb[n2]);

                    distxy =(x3 - kseedsx[n2])*(x3 - kseedsx[n2]) + (y3 - kseedsy[n2])*(y3 - kseedsy[n2]);

                    dist += distxy*invwt;
                    if(dist < distvec[index2])
                    {
                        distvec[index2] = dist;
                        klabels[index2]  = n2;
                    }
                }
            }
        }

        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        //-----------------------------------------------------------------
        for(k1 = 0; k1 < numk; k1++)
        {
            sigmal[k1] = 0;
            sigmaa[k1] = 0;
            sigmab[k1] = 0;
            sigmax[k1] = 0;
            sigmay[k1] = 0;
            clustersize[k1] = 0;
        }

        ind = 0;
        for( r = 0; r < height; r++ )
        {
            for( c = 0; c < width; c++ )
            {
                if(klabels[ind] >= 0)
                {
                    sigmal[klabels[ind]] += lvec[ind];
                    sigmaa[klabels[ind]] += avec[ind];
                    sigmab[klabels[ind]] += bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    clustersize[klabels[ind]] += 1.0;
                }
                ind++;
            }
        }

        {for( k1 = 0; k1 < numk; k1++ )
        {
            if( clustersize[k1] <= 0 ) clustersize[k1] = 1;
            inv[k1] = 1.0/clustersize[k1];//computing inverse now to multiply, than divide later
        }}

        {for( k1 = 0; k1 < numk; k1++ )
        {
                kseedsl[k1] = sigmal[k1]*inv[k1];
                kseedsa[k1] = sigmaa[k1]*inv[k1];
                kseedsb[k1] = sigmab[k1]*inv[k1];
                kseedsx[k1] = sigmax[k1]*inv[k1];
                kseedsy[k1] = sigmay[k1]*inv[k1];

        }}
    }

    //------------------------------------------------------------------------
    // Enforce connectivity
    //----------------------------------------------------------------------

    int i3,j3,k3;
    int n3,c3,count;
    int x4,y4;
    int ind3;
    int oindex, adjlabel;
    int label;
    const int dx4[4] = {-1,  0,  1,  0};
    const int dy4[4] = { 0, -1,  0,  1};
    const int SUPSZ = sz/numSuperpixels;
    int* nlabels;
    nlabels = (int*) malloc (sz*4);//corrected labels after enforcing connectivity
    int* xvec = (int*) malloc (SUPSZ*10*4);
    int* yvec = (int*) malloc (SUPSZ*10*4);

    for( i3 = 0; i3 < sz; i3++ ) nlabels[i3] = -1;
    oindex = 0;
    adjlabel = 0;//adjacent label
    label = 0;


    for( j3 = 0; j3 < height; j3++ )
    {
            for( k3 = 0; k3 < width; k3++ )
            {
                    if( 0 > nlabels[oindex] )
                    {
                            nlabels[oindex] = label;
                            //--------------------
                            // Start a new segment
                            //--------------------
                            xvec[0] = k3;
                            yvec[0] = j3;
                            //-------------------------------------------------------
                            // Quickly find an adjacent label for use later if needed
                            //-------------------------------------------------------
                            {for( n3 = 0; n3 < 4; n3++ )
                            {
                                    int x4 = xvec[0] + dx4[n3];
                                    int y4 = yvec[0] + dy4[n3];
                                    if( (x4 >= 0 && x4 < width) && (y4 >= 0 && y4 < height) )
                                    {
                                            int nindex = y4*width + x4;
                                            if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                                    }
                            }}

                            count = 1;
                            for( c3 = 0; c3 < count; c3++ )
                            {
                                    for( n3 = 0; n3 < 4; n3++ )
                                    {
                                            x4 = xvec[c3] + dx4[n3];
                                            y4 = yvec[c3] + dy4[n3];

                                            if( (x4 >= 0 && x4 < width) && (y4 >= 0 && y4 < height) )
                                            {
                                                    int nindex = y4*width + x4;

                                                    if( 0 > nlabels[nindex] && klabels[oindex] == klabels[nindex] )
                                                    {
                                                            xvec[count] = x4;
                                                            yvec[count] = y4;
                                                            nlabels[nindex] = label;
                                                            count++;
                                                    }
                                            }

                                    }
                            }
                            //-------------------------------------------------------
                            // If segment size is less then a limit, assign an
                            // adjacent label found before, and decrement label count.
                            //-------------------------------------------------------
                            if(count <= SUPSZ >> 2)
                            {
                                    for( c3 = 0; c3 < count; c3++ )
                                    {
                                            ind3 = yvec[c3]*width+xvec[c3];
                                            nlabels[ind3] = adjlabel;
                                    }
                                    label--;
                            }
                            label++;
                    }
                    oindex++;
            }
    }

    //-----------------------------------------------------------------------------------
    // Assign output labels
    //----------------------------------------------------------------------------------
    int finalNumberOfLabels = label;
    image2d <uint8_t> outlabels(D);
    std::vector< std::vector<unsigned> > listpixel(finalNumberOfLabels);
    std::vector< std::vector<V> > listcolor(finalNumberOfLabels);
    std::vector<V>  listcolor_median(finalNumberOfLabels);

    image2d <uint8_t> adjcMatrix(finalNumberOfLabels,finalNumberOfLabels);
    extension::fill (adjcMatrix, 0);
    image2d <V> Seg_Img(D);


    int k4, j4, n4;

    for(j4 = 0; j4 < width; j4++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(k4 = 0; k4 < height; k4++)
        {
            i = k4*width+j4;
            outlabels(point2d(k4,j4)) = nlabels[i];
            listpixel[nlabels[i]].push_back(i);
            listcolor[nlabels[i]].push_back(Img(point2d(k4,j4)));
            yvec[0] = k4;
            xvec[0] = j4;

            for( n4 = 0; n4 < 4; n4++ )
            {
                int x5 = xvec[0] + dx4[n4];
                int y5 = yvec[0] + dy4[n4];
                if( (x5 >= 0 && x5 < width) && (y5 >= 0 && y5 < height) )
                {
                    int nindex = y5*width + x5;
                    adjcMatrix(point2d(nlabels[i],nlabels[nindex])) = 255;
                }
            }
        }
    }


    // Find median color
    std::cout << finalNumberOfLabels << std::endl;
    for (int e = 0; e< finalNumberOfLabels; e++)
    {
        std::partial_sort(listcolor[e].begin(), listcolor[e].begin() + listcolor[e].size()/2+1, listcolor[e].end(),lexicographicalorder_less<V>());
        listcolor_median[e] = listcolor[e][int(listcolor[e].size()/2)];
    }


    for(j4 = 0; j4 < width; j4++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(k4 = 0; k4 < height; k4++)
        {
            i = k4*width+j4;
            Seg_Img(point2d(k4,j4)) = listcolor_median[nlabels[i]];
        }
    }




    //io::imsave(adjcMatrix, "vkl.ppm");
    io::imsave(Seg_Img, "dkm.ppm");

    double stop_s=clock();

    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;


    image2d<V> f = addborder(Seg_Img, lexicographicalorder_less<value_t>());
    //image2d<V> f = ima;
    io::imsave(f, "addborder.png");

    image2d<V> F = interpolate_k2(Seg_Img);




    image2d<V> ima_compo(F.domain() );







    // COmpute the color tree of shapes T
    tree_t T;
    {
      // 1. Compute the marginal ToS and filter them if necessary.
      tree_t t[NTREE];
      tbb::parallel_for(0, (int)NTREE, [&t,&Seg_Img,a0](int i){
          t[i] = morpho::cToS(imtransform(Seg_Img, [i](value_t x) { return x[i]; }), c4);
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

    stop_s=clock();

    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

          const unsigned N1 = T.nodes().size();
          std::cout << "N1  "  << N1 << std::endl;

          auto& nodes1 = T._get_data()->m_nodes;
          auto& S     = T._get_data()->m_S;
          auto& pmap  = T._get_data()->m_pmap;


          point2d startpoint= point2d(0, 0);
          std::cout <<  "start " <<pmap(startpoint)    << std::endl;

          box2d Dom = pmap.domain();
          image2d<V> distancemap(Dom );
          image2d<V> min_pixel(Dom );
          image2d<V> max_pixel(Dom );


          image2d<uint8_t> distancemap1(Dom );

          unsigned start = pmap(startpoint);

      unsigned nodenumbers = 0;


      std::cout << "node numbers   "  << T.nodes().size() << std::endl;
      nodenumbers = T.nodes().size();



      std::vector<V>  min_node(nodenumbers);
      std::vector<V>  max_node(nodenumbers);
      std::vector<V>  dist1(nodenumbers);



      std::vector< std::vector<unsigned> > children(nodenumbers);

      //int enqueue[nodenumbers] = {0};
      std::vector<int>  enqueue(nodenumbers);
      std::fill (enqueue.begin(),enqueue.end(),0);
      enqueue[0]  = 1;


      std::cout << "initiation   "  << std::endl;

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
                          dist1[x.id()][i] = max_node[x.id()][i] - min_node[x.id()][i];
                  }

                  unsigned q = x.get_parent_id();
                  children[q].push_back(x.id());
          }


          std::queue<unsigned>  Q;

      Q.push(start);
      enqueue[start] = 1;

      //auto x1 = T.get_node(start);
      //
      //std::cout << "x1.id  "  <<x1.id()  << std::endl;


          //unsigned par_p = (T.get_node(start).get_parent_id());
          //std::cout << "par_p   "  << par_p  << std::endl;

          std::cout <<  "propagation "  << std::endl;

      while(not Q.empty())
          {
                  unsigned nodeid = Q.front();
                  Q.pop();

                  for (int i = 0; i< 3; ++i)
                  {
                          dist1[nodeid][i] = max_node[nodeid][i] - min_node[nodeid][i];
                  }

                  unsigned par_p = (T.get_node(nodeid).get_parent_id());
                  if (enqueue[par_p] == 0)
                  {
                          for (int i = 0; i< 3; ++i)
                          {
                                  if (min_node[par_p][i] > min_node[nodeid][i])
                                          min_node[par_p][i] = min_node[nodeid][i];
                                  if (max_node[par_p][i] < max_node[nodeid][i])
                                          max_node[par_p][i] = max_node[nodeid][i];
                          }
                          Q.push(par_p);
                          enqueue[par_p] =1;
                  }
                  const unsigned nchildren = children[nodeid].size();

                  for (unsigned j= 0; j< nchildren ; ++j)
                  {
                          unsigned child = children[nodeid][j];
                          if (enqueue[child] == 0)
                          {
                                  for (int i = 0; i< 3; ++i)
                                  {
                                          if (min_node[child][i] > min_node[nodeid][i])
                                                  min_node[child][i] = min_node[nodeid][i];
                                          if (max_node[child][i] < max_node[nodeid][i])
                                                  max_node[child][i] = max_node[nodeid][i];
                                  }
                                  Q.push(child);
                                  enqueue[child] =1;
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

                  }

          }

          stop_s=clock();

          std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

      io::imsave(distancemap, "dmap.png");
      io::imsave(distancemap1, "dmap_gray.png");



}
