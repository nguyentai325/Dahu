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

#include <utility>
#include <algorithm>
#include <vector>

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"


#include "range_sp.hpp"








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

    bool is_border (point2d p, unsigned height , unsigned width)
    {
        if (p[0] == 0 or p[0]== height -1 or p[1] == 0 or p[1] == width)
            return true;
        else
            return false;

    }



    uint8 upper_level_next_to_lcur(std::vector<std::queue<unsigned> > & q, uint8 l_cur, bool& found)
    {
        uint8 v = l_cur;
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


    uint8 lower_level_next_to_lcur(std::vector<std::queue<unsigned> > & q, uint8 l_cur, bool& found)
    {
        uint8 v = l_cur;
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

    uint8 level_next_to_lcur(std::vector<std::queue<unsigned> > & q, uint8 l_cur)
    {
        uint8 l_;
        bool found;

        bool up = int(2. * std::rand() / (RAND_MAX + 1.));
        if (up)
        {
            //std::cout << "up  "  << std::endl;
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
            //std::cout << "down  "  << std::endl;

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

    unsigned priority_pop(std::vector<std::queue<unsigned> > &q, uint8& l_cur, int &n_element_in_queue )
    // modify q, and sometimes l_cur
    {
        if (q[l_cur].empty())
        {
            uint8 l_ = level_next_to_lcur(q, l_cur);  // such as q[l_] is not empty
            if (q[l_].empty())
            std::abort();
            l_cur = l_;
        }

        unsigned sp = q[l_cur].front();
        q[l_cur].pop();
        n_element_in_queue = n_element_in_queue -1;
        return sp;
    }

    uint8 priority_push(std::vector<std::queue<unsigned> >& q, unsigned sp,
          std::vector< range_sp<uint8> > & U,
          uint8 l_cur, int &n_element_in_queue)
    // modify q
    {
        uint8
        lower = U[sp].lower,
        upper = U[sp].upper,
        l_;
        if (lower > l_cur)
            l_ = lower;
        else if (upper < l_cur)
            l_ = upper;
        else
            l_ = l_cur;

        q[l_].push(sp);
        n_element_in_queue = n_element_in_queue + 1;

        return l_;
    }

    unsigned find_root(std::vector<unsigned> &zpar, unsigned x)
    // modify zpar
    {
      if (zpar[x] == x)
        return x;
      else
      {
        zpar[x] = find_root(zpar, zpar[x]);
        return zpar[x];
      }
    }


    void do_union(unsigned p_, unsigned r_,
          std::vector<unsigned> &zpar,
          std::vector<unsigned> &rank,
          std::vector<unsigned> &last)
    // modify zpar, rank, and last
    {
        if (rank[p_] > rank[r_])
        {
          // new root is p_
          zpar[r_]= p_;
          if (last[r_] < last[p_])
            last[p_] = last[r_];
        }
          else
        {
          // new root is r_
          zpar[p_] = r_;
          if (last[p_] < last[r_])
            last[r_] = last[p_];
          if (rank[p_] == rank[r_])
            rank[r_] = rank[r_] + 1;
        }
    }



    struct tree
    {
        std::vector<unsigned> S;
        std::vector<unsigned> parent;
        std::vector<uint8> Ub;

//        tree (const unsigned finalNumberOfLabels1)
//        {
//            std::vector<unsigned> S;
//            parent = std::vector<unsigned> (finalNumberOfLabels1);
//            Ub = std::vector<uint8> (finalNumberOfLabels1);
////            std::vector<unsigned> parent;
////            std::vector<uint8> Ub ;

//        }

        bool is_root(unsigned sp) const
        {
            return parent[sp] == sp;
        }

        bool is_representative(unsigned sp) const
        {
            return is_root(sp) || Ub[parent[sp]] != Ub[sp];
        }

        unsigned get_representative(const unsigned sp) const
        {
          unsigned x = sp;
          while (! is_representative(x))
            x = parent[x];
          return x;
        }

//        unsigned get_representative(const unsigned& p) const
//        {
//          unsigned x = p;
//          while (! is_representative(x))
//            x = parent[x];
//          return x;
//        }

        void canonicalize()
        {
            const unsigned N = S.size();
            for (unsigned i = 0; i < N; ++i)
            {
                unsigned sp = S[i];  // p goes from root to leaves
                unsigned q = parent[sp];

                // q has always been processed so parent(q) is a
                // representative node:
                if (! is_representative(parent[q]))
                    std::abort();

                // if p does not point towards a representative node then
                // skip the node above p (i.e., p points to q's parent,
                // which is a representative node)

                if (Ub[parent[q]] == Ub[q])
                    parent[sp] = parent[q];
            }
        }

        std::vector <unsigned>  representative_node (tree T)
        {
            std::vector<unsigned> representation;
            for (unsigned i = 0; i < T.S.size(); ++i)
            {

                if (T.is_representative(T.S[i]))
                {
                    representation.push_back(T.S[i]);
                }
            }
            return representation;
        }

        std::vector<unsigned> back_propagate(tree T, std::vector<unsigned> A)
        {
            for (unsigned i = 0; i < T.S.size(); ++i) // from root to leaves
            {
                if (! T.is_representative(T.S[i]))
                    A[T.S[i]] = A[T.parent[T.S[i]]];
            }
        }





        // back-propagate


        int least_common_ancestor(int A, int B, std::vector <unsigned >  depth ,tree T )
        {
            if (A == -1)
                return B;

            if (B == -1)
                return A;

            while(depth[A] < depth[B])
                B = T.parent[B];
            while(depth[B] < depth[A])
                A = T.parent[A];
            while(A!= B)
            {
                A = T.parent[A];
                B = T.parent[B];
            }
            return A;
        }





    };


    int node_to_index(std::vector<unsigned>  node_graph, unsigned node)
    {
        for(int i =0; i < node_graph.size() ; i++)
        {
            if (node_graph[i]  == node)
                return i;
        }
        return 0;
    }


}




void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input[rgb]  num \n"
    "α₀\tGrain filter size before merging trees (0 to disable)\n"
    "α₁\tGrain filter size on the color ToS (0 to disable)\n"
    "λ\tMumford-shah regularisation weight (e.g. 5000)\n";
  std::exit(1);
}



int main(int argc, char** argv)
{
    if (argc < 2)
    usage(argv);

    double start_s=clock();

    const char* input_path = argv[1];
    int numSuperpixels = std::atoi(argv[2]);
    std::string fileName = input_path;

    using namespace boost;

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

    std::cout << "height"  <<height<< std::endl;
    std::cout << "width"  << width << std::endl;

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

    const int dx8[8] = {-1,  0,  1,  0, -1, 1,  1, -1};
    const int dy8[8] = { 0, -1,  0,  1, -1, 1, -1,  1};

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
    std::cout << "number of super pixel"  << finalNumberOfLabels  << std::endl;
    image2d <uint8_t> outlabels(D);
    std::vector< std::vector<unsigned> > listpixel(finalNumberOfLabels+1);
    std::vector< std::vector<V> > listcolor(finalNumberOfLabels+1);
    std::vector<V>  listcolor_median(finalNumberOfLabels+1);
    std::vector<uint8_t>  listcolor_median_gray(finalNumberOfLabels+1);

    std::vector<uint8_t>  listcolor_median_red(finalNumberOfLabels+1);
    std::vector<uint8_t>  listcolor_median_green(finalNumberOfLabels+1);
    std::vector<uint8_t>  listcolor_median_blue(finalNumberOfLabels+1);





    std::vector<std::vector<unsigned> > neighbor_sp(finalNumberOfLabels+1);
    std::vector<unsigned> sp_border;

    image2d <uint8_t> adjcMatrix(finalNumberOfLabels+1,finalNumberOfLabels+1);
    extension::fill (adjcMatrix, 0);
    image2d <V> Seg_Img(D);
    image2d <uint8_t>  Seg_Img_gray(D);


    int k4, j4, n4;

    for(j4 = 0; j4 < height; j4++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(k4 = 0; k4 < width; k4++)
        {
            i = j4*width+k4;
            outlabels(point2d(j4,k4)) = nlabels[i];
            listpixel[nlabels[i]].push_back(i);
            listcolor[nlabels[i]].push_back(Img(point2d(j4,k4)));
            //yvec[0] = j4;
            //xvec[0] = k4;

            if (is_border(point2d(j4,k4),height,width))
            {
                 sp_border.push_back(nlabels[i]);
            }

            // CONNECTION C8
            for( n4 = 0; n4 < 8; n4++ )
            {
                int x5 = k4 + dx8[n4];
                int y5 = j4 + dy8[n4];
                if( (x5 >= 0 && x5 < width) && (y5 >= 0 && y5 < height) )
                {
                    int nindex = y5*width + x5;
                    adjcMatrix(point2d(nlabels[i],nlabels[nindex])) = 255;
                    if (nlabels[i] != nlabels[nindex])
                    {
                        if (neighbor_sp[nlabels[i]].empty())
                            neighbor_sp[nlabels[i]].push_back(nlabels[nindex]);
                        else
                        {
                            bool found = false;
                            for (int t = 0; t < neighbor_sp[nlabels[i]].size() ; t++)
                            {
                                if (nlabels[nindex] == neighbor_sp[nlabels[i]][t])
                                {
                                    found = true;
                                    break;
                                }
                            }
                            if (found == false)
                                neighbor_sp[nlabels[i]].push_back(nlabels[nindex]);
                        }
                    }
                }
            }
        }
    }







    // Find median color
    for (int e = 0; e< finalNumberOfLabels; e++)
    {
        std::partial_sort(listcolor[e].begin(), listcolor[e].begin() + listcolor[e].size()/2+1, listcolor[e].end(),lexicographicalorder_less<V>());
        listcolor_median[e] = listcolor[e][int(listcolor[e].size()/2)];
        listcolor_median_gray[e] = 0.2989 * listcolor_median[e][0] + 0.5870 * listcolor_median[e][1] + 0.1140 * listcolor_median[e][2];
        listcolor_median_red[e] =  listcolor_median[e][0] ;
        listcolor_median_green[e] = listcolor_median[e][1];
        listcolor_median_blue[e] = listcolor_median[e][2];


    }



    // superpixel on the border of the image.

    std::vector<unsigned>::iterator ip;

    // Sorting the array
    std::sort(sp_border.begin(), sp_border.end());

    // Using std::unique
    ip = std::unique(sp_border.begin(), sp_border.end());

    // Resizing the vector so as to remove the undefined terms
    sp_border.resize(std::distance(sp_border.begin(), ip));

    std::vector<V> listcolor_border;
    uint8 listcolor_border_gray;
    int label_border = finalNumberOfLabels;
    std::cout << "label border  "  << label_border << std::endl;


    for (int t = 0; t < sp_border.size() ; t++)
    {
        listcolor_border.push_back(listcolor_median[sp_border[t]]);
        adjcMatrix(point2d(sp_border[t],finalNumberOfLabels)) = 255;
        adjcMatrix(point2d(finalNumberOfLabels,sp_border[t])) = 255;
        neighbor_sp[sp_border[t]].push_back(finalNumberOfLabels);
        neighbor_sp[finalNumberOfLabels].push_back(sp_border[t]);
    }

    std::partial_sort(listcolor_border.begin(), listcolor_border.begin() + listcolor_border.size()/2+1, listcolor_border.end(),lexicographicalorder_less<V>());
    listcolor_median[finalNumberOfLabels] = listcolor_border[int(listcolor_border.size()/2)];
    std::cout << int(listcolor_median[finalNumberOfLabels][0])  << int(listcolor_median[finalNumberOfLabels][1])  << int(listcolor_median[finalNumberOfLabels][2])  << std::endl;
    listcolor_median_gray[finalNumberOfLabels] = 0.2989 * listcolor_median[finalNumberOfLabels][0] + 0.5870 * listcolor_median[finalNumberOfLabels][1] + 0.1140 * listcolor_median[finalNumberOfLabels][2];
    listcolor_median_red[finalNumberOfLabels]   =  listcolor_median[finalNumberOfLabels][0] ;
    listcolor_median_green[finalNumberOfLabels] =  listcolor_median[finalNumberOfLabels][1] ;
    listcolor_median_blue[finalNumberOfLabels]  =  listcolor_median[finalNumberOfLabels][2] ;


    finalNumberOfLabels = finalNumberOfLabels + 1;



    for(j4 = 0; j4 < height; j4++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(k4 = 0; k4 < width; k4++)
        {
            i = j4*width+k4;
            Seg_Img(point2d(j4,k4)) = listcolor_median[nlabels[i]];
            Seg_Img_gray(point2d(j4,k4)) = listcolor_median_gray[nlabels[i]];

        }
    }

    io::imsave(Seg_Img, "seg.ppm");






    ////////////////////////////////////////////////  Tree of shapes  //////////////////////////////////////////////


    // graph of super-pixel

    //create an -undirected- graph type, using vectors as the underlying containers
    //and an adjacency_list as the basic representation
    typedef adjacency_list<vecS, vecS, undirectedS> UndirectedGraph;

    //An edge is just a connection between two vertitices. Our verticies above
    //are an enum, and are just used as integers, so our edges just become
    //a std::pair<int, int>
    typedef std::pair<int, int> Edge;
    UndirectedGraph g(finalNumberOfLabels);

    //std::vector<Edge> edgeVec;
    for(int t = 0; t < neighbor_sp.size() ; t++)
    {

        std::sort(neighbor_sp[t].begin(), neighbor_sp[t].end());
        for(int u = 0; u < neighbor_sp[t].size() ; u++ )
        {
            if (neighbor_sp[t][u] > t)
                add_edge(t, neighbor_sp[t][u], g);

        }
    }





    typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;
    IndexMap index_graph1 = get(vertex_index, g);
    typename graph_traits<UndirectedGraph>::adjacency_iterator ai1;
    typename graph_traits<UndirectedGraph>::adjacency_iterator ai_end1;


    typedef graph_traits<UndirectedGraph>::vertex_descriptor Vertex;

    // get the property map for vertex indices


    //std::cout << "vertices(g) = ";
    typedef graph_traits<UndirectedGraph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;



    //   find triangle nodes in the graph

    //std::cout << "edges(g) = ";

    typedef std::array<unsigned, 3> triangle;

    std::vector<triangle> list_triangle;
    graph_traits<UndirectedGraph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
//        std::cout << "(" << index_graph1[source(*ei, g)]
//                  << "," << index_graph1[target(*ei, g)] << ") ";

        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
        {
          Vertex v = *vp.first;
          triangle tri_check;
          tri_check[0] = index_graph1[source(*ei, g)];
          tri_check[1] = index_graph1[target(*ei, g)];
          tri_check[2] = index_graph1[v];
          std::sort (tri_check.begin(), tri_check.begin()+3);

          //std::cout << index_graph1[v] <<  " ";
          if (adjcMatrix(point2d(index_graph1[source(*ei, g)],v)) == 255 and adjcMatrix(point2d(index_graph1[target(*ei, g)],v)) == 255 and v != index_graph1[source(*ei, g)] and v != index_graph1[target(*ei, g)])
          {
              unsigned check_tong = 0;
              for (int t = 0; t < list_triangle.size() ; t++)
              {
                    triangle tri = list_triangle[t];
                    std::sort (tri.begin(), tri.begin()+3);
                    unsigned check = 0;
                    for (int s = 0; s < 3; s++)
                    {
                        if (tri[s] == tri_check[s])
                            check = check +1;
                    }
                    if (check < 3)
                        check_tong = check_tong +1;


              }
              if (check_tong == list_triangle.size())
                    list_triangle.push_back(tri_check);
          }

        }


    }

    //std::cout << "triangle number "  << list_triangle.size()  << std::endl;


    //Now we can initialize our graph using iterators from our above vector
    //UndirectedGraph g(edgeVec.begin(), edgeVec.end(), N);


    //Ok, we want to see that all our edges are now contained in the graph
    typedef graph_traits<UndirectedGraph>::edge_iterator edge_iterator;






//    std::vector<tree> T(3);
    tree T[3];

    unsigned finalNumberOfLabels1 =0;
    for (int channel = 0; channel < 3 ; channel++)
    {
        //// Create inter-superpixel (immerse)
        label = finalNumberOfLabels;

        std::cout << "channel  "  << channel << std::endl;

        typedef range_sp<uint8> R;
        std::vector< range_sp<uint8> >  range_color_sp(finalNumberOfLabels + num_edges(g) );

        for(int t = 0; t <finalNumberOfLabels ; t++)
        {
            if (channel == 0)
            {
                range_color_sp[t] = listcolor_median_red[t];

                std::sort(neighbor_sp[t].begin(), neighbor_sp[t].end());

                for(int u = 0; u < neighbor_sp[t].size() ; u++ )
                {
                    if (neighbor_sp[t][u] > t)
                    {
                        add_edge(t, label, g);
                        add_edge(neighbor_sp[t][u], label, g);
                        remove_edge(t,neighbor_sp[t][u],g);
                        range_color_sp[label] = R{std::min(listcolor_median_red[t],listcolor_median_red[neighbor_sp[t][u]]), std::max(listcolor_median_red[t],listcolor_median_red[neighbor_sp[t][u]])};
                        //std::cout << int (range_color_sp[label].lower) << "   " << int(range_color_sp[label].upper) << std::endl;
                        label = label +1;
                    }
                }
            }
            else if (channel == 1)
            {
                range_color_sp[t] = listcolor_median_green[t];

                std::sort(neighbor_sp[t].begin(), neighbor_sp[t].end());

                for(int u = 0; u < neighbor_sp[t].size() ; u++ )
                {

                    if (neighbor_sp[t][u] > t)
                    {
                        add_edge(t, label, g);
                        add_edge(neighbor_sp[t][u], label, g);
                        remove_edge(t,neighbor_sp[t][u],g);
                        range_color_sp[label] = R{std::min(listcolor_median_green[t],listcolor_median_green[neighbor_sp[t][u]]), std::max(listcolor_median_green[t],listcolor_median_green[neighbor_sp[t][u]])};
                        //std::cout << int (range_color_sp[label].lower) << "   " << int(range_color_sp[label].upper) << std::endl;
                        label = label +1;
                    }


                }
            }
            else
            {
                range_color_sp[t] = listcolor_median_blue[t];

                std::sort(neighbor_sp[t].begin(), neighbor_sp[t].end());

                for(int u = 0; u < neighbor_sp[t].size() ; u++ )
                {
                    if (neighbor_sp[t][u] > t)
                    {
                        add_edge(t, label, g);
                        add_edge(neighbor_sp[t][u], label, g);
                        remove_edge(t,neighbor_sp[t][u],g);
                        range_color_sp[label] = R{std::min(listcolor_median_blue[t],listcolor_median_blue[neighbor_sp[t][u]]), std::max(listcolor_median_blue[t],listcolor_median_blue[neighbor_sp[t][u]])};
                        //std::cout << int (range_color_sp[label].lower) << "   " << int(range_color_sp[label].upper) << std::endl;
                        label = label +1;
                    }
                }
            }
        }



        finalNumberOfLabels1 = label;

//        std::cout << finalNumberOfLabels1 << std::endl;
//        std::cout << finalNumberOfLabels << std::endl;

        typedef property_map<UndirectedGraph, vertex_index_t>::type IndexMap;
        IndexMap index_graph = get(vertex_index, g);
        typename graph_traits<UndirectedGraph>::adjacency_iterator ai;
        typename graph_traits<UndirectedGraph>::adjacency_iterator ai_end;





        for (int t = 0; t < list_triangle.size() ; t++)
        {
            triangle tri = list_triangle[t];

            std::vector <unsigned > adja_1;
            for (boost::tie(ai, ai_end) = adjacent_vertices(tri[0], g); ai != ai_end; ++ai)
                 adja_1.push_back(index_graph[*ai]);

            std::vector <unsigned > adja_2;
            for (boost::tie(ai, ai_end) = adjacent_vertices(tri[1], g); ai != ai_end; ++ai)
                 adja_2.push_back(index_graph[*ai]);

            std::vector <unsigned > adja_3;
            for (boost::tie(ai, ai_end) = adjacent_vertices(tri[2], g); ai != ai_end; ++ai)
                 adja_3.push_back(index_graph[*ai]);




            std::sort (adja_1.begin(), adja_1.begin() + adja_1.size());
            std::sort (adja_2.begin(), adja_2.begin() + adja_2.size());
            std::sort (adja_3.begin(), adja_3.begin() + adja_3.size());

            std::vector<unsigned> common(adja_1.size() + adja_2.size());
            std::vector<unsigned>::iterator it;
            it=std::set_intersection (adja_1.begin(), adja_1.begin() +adja_1.size(), adja_2.begin(),  adja_2.begin() + adja_2.size(), common.begin());
            common.resize(it-common.begin());

            std::vector<unsigned> common1(adja_1.size() + adja_3.size());
            std::vector<unsigned>::iterator it1;
            it1=std::set_intersection (adja_1.begin(), adja_1.begin() +adja_1.size(), adja_3.begin(),  adja_3.begin() + adja_3.size(), common1.begin());
            common1.resize(it1-common1.begin());

            std::vector<unsigned> common2(adja_2.size() + adja_3.size());
            std::vector<unsigned>::iterator it2;
            it2=std::set_intersection (adja_2.begin(), adja_2.begin() +adja_2.size(), adja_3.begin(),  adja_3.begin() + adja_3.size(), common2.begin());
            common2.resize(it2-common2.begin());


            if (common[0] < common1[0])
                add_edge(common[0], common1[0], g);
            if (common[0] < common2[0])
                add_edge(common[0], common2[0], g);
            if (common1[0] < common2[0])
                add_edge(common1[0], common2[0], g);

        }


//        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
//        {
//            std::cout << "(" << index_graph1[source(*ei, g)]
//                      << "," << index_graph1[target(*ei, g)] << ") ";
//        }






        //// sort super pixel
        ///
        /// tree
        ///
        ///


        T[channel].parent = std::vector<unsigned> (finalNumberOfLabels1);
        T[channel].Ub = std::vector<uint8> (finalNumberOfLabels1);


        // q and deja_vu

        std::vector<std::queue<unsigned> > q(256);

        std::vector<bool>  deja_vu(finalNumberOfLabels1);
        std::fill (deja_vu.begin(),deja_vu.end(),false);
        std::vector<uint8> U_sp(finalNumberOfLabels1);

        // distance map
        typedef std::pair<uint8,uint8> pair_t;
        std::vector<pair_t> mm(finalNumberOfLabels1);


        // auxiliary data, just for debug purpose

        std::vector<bool> done(finalNumberOfLabels1);
        std::fill (done.begin(),done.end(),false);


        // initialization
        unsigned sp_oo = label_border;
        uint8 l_cur = range_color_sp[sp_oo].lower;
        std::vector< unsigned > sort_sp;
        int n_element_in_queue = 0;

        q[l_cur].push(sp_oo);
        n_element_in_queue = n_element_in_queue + 1;
        deja_vu[sp_oo] = true;





        mm[sp_oo] = pair_t(l_cur, l_cur);

        while (n_element_in_queue != 0)
        {
            unsigned sp = priority_pop(q, l_cur, n_element_in_queue);
            sort_sp.push_back(sp);
            //std::cout << "sp  "  << sp   << "  "<< int(l_cur)<< std::endl;



            T[channel].Ub[sp] = l_cur;
            //std::cout << int(l_cur) << std::endl;
            T[channel].S.push_back(sp);


            for (boost::tie(ai, ai_end) = adjacent_vertices(sp, g); ai != ai_end; ++ai)
            {

                unsigned n_sp = index_graph[*ai];

                if (deja_vu[n_sp] == false)
                {
                    //std::cout << "n_sp  "  << n_sp << std::endl;

                    uint8 l_ = priority_push(q, n_sp, range_color_sp, l_cur, n_element_in_queue);
                   // std::cout << int(range_color_sp[n_sp].lower) << "  "  << int(range_color_sp[n_sp].upper) << " "<< int(l_)  <<"  "  << int(l_cur) << std::endl;
                    U_sp[n_sp] = l_;
                    deja_vu[n_sp] = true;
                    mm[n_sp] = mm[sp];
                    if (l_ < mm[n_sp].first)
                        mm[n_sp].first = l_;
                    if (l_ > mm[n_sp].second)
                        mm[n_sp].second = l_;
                }
            }
            done[sp] = true;
        }

        for (int t = 0 ; t< finalNumberOfLabels1 ; t++)
        {
            if (done[t] == false)
            {
                std::cout << "sai vkl  "  << std::endl;
            }
        }


        // compute dahu map

        image2d<uint8_t>  dmap(D);



        for (int t = 0; t< finalNumberOfLabels; t++)
        {
            for(int e = 0; e < listpixel[t].size() ; e++)
            {
                unsigned py = listpixel[t][e]/width;
                unsigned px = listpixel[t][e]%width;
                dmap(point2d(py,px)) = mm[t].second - mm[t].first;

            }
        }


        io::imsave(dmap, "vkl.ppm");




        // union find


        // shortcuts



        // auxiliary data

        unsigned undef = INT_MAX;

        std::vector<unsigned> zpar(finalNumberOfLabels1);
        std::fill (zpar.begin(),zpar.end(),undef);

        std::vector<unsigned> rank(finalNumberOfLabels1), last(finalNumberOfLabels1);
        std::fill (rank.begin(),rank.end(),0);
        std::fill (last.begin(),last.end(),0);

    //    // end of auxiliary data



        // only used for display:
        std::vector<bool> vu(finalNumberOfLabels1);
        std::fill (vu.begin(),vu.end(),false);

        std::cout <<"size tree  " <<T[channel].S.size()  << std::endl;

        for (int t = T[channel].S.size() -1  ; t >= 0; --t)
        {
            unsigned sp = T[channel].S[t];
            //std::cout << "sp  "  << sp << std::endl;
            T[channel].parent[sp] = sp;
            zpar[sp] = sp;
            last[sp] = t;

            for (boost::tie(ai, ai_end) = adjacent_vertices(sp, g); ai != ai_end; ++ai)
            {
                unsigned n_sp = index_graph[*ai];
                if (zpar[n_sp] != undef)
                {
                    //std::cout << "n_sp  "<<n_sp << std::endl;

                    unsigned p_ = find_root(zpar,sp);
                    unsigned r_ = find_root(zpar,n_sp);
                    if (r_ != p_)
                    {
                        unsigned r = T[channel].S[last[r_]];
                        T[channel].parent[r] = sp;
                        //std::cout << r  << "  has parent  "  << sp << std::endl;
                        do_union(p_,r_,zpar,rank,last);
                        //zpar[r_] = sp;
                    }
                }
            }
            vu[sp] = true;

        }


        // canonicalize

        T[channel].canonicalize();
    }

    ///////////////////////////////  end tree of shape  /////////////////////////////////////////////






    /////////////////////////////////////////////  GRAPH PRECOMPUTATION //////////////////////////////////////////

    ////// check if representative node is 2faces or not
    //// if node we find the first 2D children and swap it with the representative node ...
    ///


    for(int channel = 0 ; channel < 3 ; channel ++)
    {
        std::vector< std::vector<unsigned> > children(T[channel].S.size());

        std::vector< std::vector<unsigned> > index_children(T[channel].S.size());


        //std::cout << "tree  " << channel  << std::endl;
        for (unsigned i = 0; i < T[channel].S.size(); i++)
        {
            unsigned p = T[channel].S[i];
            unsigned q = T[channel].parent[p];
            children[q].push_back(p);
            index_children[q].push_back(i);
//            if (T[channel].is_representative(p))
//                std::cout << p  << "  "<< q << std::endl;

        }






        for (unsigned i = 0; i < T[channel].S.size(); i++)
        {
            unsigned p = T[channel].S[i];
            if (T[channel].is_representative(p)  and p  >= finalNumberOfLabels)   // not 2D faces
            {
                for(unsigned j = 0; j < children[p].size(); j++)
                {
                    unsigned child = children[T[channel].S[i]][j];

                    //std::cout << p   << " is parent of  "  << child   << std::endl;

                    if (!T[channel].is_representative(child) and child < finalNumberOfLabels) // first 2faces
                    {
                        //std::cout << "parent  " << p <<"   child   "  << child << std::endl;
                        T[channel].S[i] = child;                      // swap child and its parent in order of S
                        T[channel].S[index_children[p][j]] = p;


                        unsigned temp_parent = T[channel].parent[p];   // swap parent relation
                        T[channel].parent[p] = child;
                        T[channel].parent[child] = temp_parent;
                        children[child].push_back(p);             // add new child to the new parent



                        for (unsigned k = 0; k < children[p].size(); k++)
                        {
                            unsigned child_child = children[p][k];

                            if (child_child != child)
                            {
                                T[channel].parent[child_child] = child;
                                children[child].push_back(child_child);
                            }

                        }
                        children[p].erase(children[p].begin(),children[p].begin()+children[p].size());

                        break;

                    }

                }
            }
        }


        //// test plot tree
//        std::cout << "channel   " << channel << std::endl;

//        for (unsigned i = 0; i < T[channel].S.size(); i++)
//        {
//            unsigned p = T[channel].S[i];
//            if (T[channel].is_representative(p) )
//            {
//                for(unsigned j = 0; j < children[p].size(); j++)
//                {
//                    unsigned child = children[T[channel].S[i]][j];

//                    std::cout << p   << " is parent of  "  << child   << std::endl;

//                }
//            }
//        }


    }




    for(int channel = 0 ; channel < 3 ; channel ++)
    {
        std::vector<unsigned> representation = T[channel].representative_node(T[channel]);
        for (int i = 0 ; i < representation.size() ; i++)
        {
            std::vector<uint8> node(T[channel].S.size());
            std::fill(node.begin(),node.end(),0);
            if (T[channel].is_representative(representation[i]))
            {
                node[representation[i]] = 255;
                for(int j= 0 ; j < T[channel].S.size(); j++)
                {
                    unsigned p = T[channel].S[j];
                    //std::cout << p << "   "  << T[0].parent[T[0].S[j]]  << std::endl;
                    if (node[p] == 0 and node[T[channel].parent[p]] == 255)
                    {
                        node[p] = node[T[channel].parent[p]];
                    }
                }

                image2d<uint8_t>  node_ima(D);

                for (int t = 0; t< finalNumberOfLabels; t++)
                {
                    for(int e = 0; e < listpixel[t].size() ; e++)
                    {
                        unsigned py = listpixel[t][e]/width;
                        unsigned px = listpixel[t][e]%width;
                        node_ima(point2d(py,px)) = node[t];
                    }
                }

                std::string a = "node_" + std::to_string(channel) + "_"+ std::to_string(representation[i]) + ".png";
                io::imsave(node_ima, a.c_str());

            }
        }
    }






    //// compute depth node for 3 trees

    unsigned depth[3][T[0].S.size()];
    memset( depth, 0, 3*T[0].S.size()*sizeof(unsigned) );
    for (unsigned channel = 0; channel < 3 ; channel++)
    {
        //std::cout << "channel  "  << channel << std::endl;
        for (unsigned i = 0; i < T[channel].S.size(); i++)
        {
            if (i == 0)
                depth[channel][T[channel].S[i]] = 0;
            else
                depth[channel][T[channel].S[i]] = depth[channel][T[channel].parent[T[channel].S[i]]] + 1;
            //std::cout << T[channel].S[i] << "   " << depth[channel][T[channel].S[i]] << std::endl;
        }
    }


    //// compute smallest closing shape



    int  SES[6][T[0].S.size()];
    memset( SES, -1, 6*T[0].S.size()*sizeof(unsigned) );

    std::cout << SES[0][0]  << std::endl;


    unsigned index_SES = 0;
    for (unsigned channel1 = 0 ; channel1 < 3 ; channel1++)
    {
        for (unsigned channel2 = 0 ; channel2 < 3 ; channel2++)
        {
            if (channel1 != channel2)
            {
                for (int i = 0; i < T[channel1].S.size(); i++)
                {
                    unsigned x = T[channel1].S[i];
                    unsigned A = T[channel1].get_representative(x);
                    unsigned X = T[channel2].get_representative(x);
                    std::vector<unsigned> extract_depth(depth[channel2], depth[channel2] + T[channel2].S.size());
                    SES[index_SES][A]  = T[channel2].least_common_ancestor(SES[index_SES][A], X, extract_depth ,T[channel2] );

                    //std::cout << T[0].S[i]  << "   " << A << "   "  << X   << "   " << SES[0][A]  << std::endl;

                }

                for (int i = T[channel1].S.size()-1; i >= 0 ; i--)
                {
                    unsigned A = T[channel1].get_representative(T[channel1].S[i]);
                    unsigned Q = T[channel1].parent[A];
                    std::vector<unsigned> extract_depth(depth[channel2], depth[channel2] + T[channel2].S.size());
                    SES[index_SES][Q] = T[channel2].least_common_ancestor(SES[index_SES][Q], SES[index_SES][A], extract_depth ,T[channel2] );
                }
                index_SES = index_SES + 1;
            }
        }
    }




//    for (int i = 0; i < T[1].S.size() ; i++)
//    {
//        std::cout << T[1].get_representative(T[1].S[i])  << "  SES  "  <<   SES[2][T[1].get_representative(T[1].S[i])]  << std::endl;
//    }


    ////////////////////////////////////////////////  GRAPH OF SHAPE COMPUTATION  /////////////////////////////////////////
    //// From 3 trees separated we construct a graph that combines 3 trees together
    ///
    ///
    // Insert the node (not twice)

    // add all nodes from the first tree, find the relation between the nodes on the old tree with the nodes on the new graph

    std::vector<unsigned> representation_T0 = T[0].representative_node(T[0]);
    std::vector<unsigned> representation_T1 = T[1].representative_node(T[1]);
    std::vector<unsigned> representation_T2 = T[1].representative_node(T[2]);

    std::vector<bool> new_node_T1(representation_T1.size());
    std::fill(new_node_T1.begin(),new_node_T1.end(), false);
    std::vector<bool> new_node_T2(representation_T2.size());
    std::fill(new_node_T2.begin(),new_node_T2.end(), false);




    // first tree

    unsigned max_size = std::max(representation_T0.size(),std::max(representation_T1.size(),representation_T2.size())) + 1;
    //std::cout << representation_T0.size()  << "  " << representation_T1.size()  << "   "  << representation_T2.size()  << std::endl;
    std::vector<unsigned> node_graph(representation_T0);
    std::vector<unsigned> node_T0_graph(max_size);
    std::vector<unsigned> node_T1_graph(max_size);
    std::vector<unsigned> node_T2_graph(max_size);
    unsigned max_graph = 0;

    for (int i = 0; i < representation_T0.size() ; i++)
    {
        //std::cout << representation_T0[i] << std::endl;
        node_T0_graph[representation_T0[i]] = representation_T0[i];
        if (max_graph < representation_T0[i])
            max_graph = representation_T0[i];
    }

    // second tree



    for (int i = 0; i < representation_T1.size() ; i++)
    {
        //std::cout << representation_T1[i]  << "   "  << SES[0][SES[2][representation_T1[i]]]  << std::endl;
        if (representation_T1[i] != SES[0][SES[2][representation_T1[i]]])
        {
            node_graph.push_back(representation_T1[i]+T[0].S.size());
            new_node_T1[representation_T1[i]] = true;
            node_T1_graph[representation_T1[i]] = representation_T1[i]+T[0].S.size();
            if (max_graph < representation_T1[i]+T[0].S.size())
                max_graph = representation_T1[i]+T[0].S.size();
        }
        else
            node_T1_graph[representation_T1[i]] = SES[2][representation_T1[i]];
    }

   // third tree



    for (int i = 0; i < representation_T2.size() ; i++)
    {
        if (representation_T2[i] == SES[1][SES[4][representation_T2[i]]])
            node_T2_graph[representation_T2[i]] = SES[4][representation_T2[i]];
        else if (representation_T2[i] == SES[3][SES[5][representation_T2[i]]])
            node_T2_graph[representation_T2[i]] = node_T1_graph[SES[5][representation_T2[i]]];
        else // (representation_T2[i] != SES[1][SES[4][representation_T2[i]]] and  representation_T2[i] != SES[3][SES[5][representation_T2[i]]])
        {
            node_graph.push_back(representation_T2[i]+2 * T[0].S.size());
            new_node_T2[representation_T2[i]] = true;
            node_T2_graph[representation_T2[i]] = representation_T2[i]+2* T[0].S.size();
            if (max_graph < representation_T2[i]+2 * T[0].S.size())
                max_graph = representation_T2[i]+2 * T[0].S.size();

        }
    }

//    // Add edges between new nodes on the graph


    UndirectedGraph graph_of_shape(node_graph.size());
    unsigned  check_edge[max_graph+1][max_graph+1];
    memset( check_edge, 0, (max_graph+1)*(max_graph+1)*sizeof(unsigned) );



    std::vector<std::vector<unsigned> > parent_graph(max_graph+1);
    std::vector<unsigned> depth_graph(max_graph+1);
    std::fill(depth_graph.begin(), depth_graph.end(), 0);




//    // first tree

    for (int i = 0; i < representation_T0.size() ; i++)
    {

        if (!T[0].is_root(representation_T0[i]))
        {
            unsigned node_t = representation_T0[i];
            add_edge(node_t, T[0].parent[node_t], graph_of_shape);
            check_edge[node_t][T[0].parent[node_t]] =1;
            parent_graph[node_t].push_back(T[0].parent[node_t]);
            depth_graph[node_t] = depth_graph[T[0].parent[node_t]] + 1;

            unsigned SES_T1 = SES[0][node_t];
            unsigned SES_T2 = SES[1][node_t];


            if (node_T1_graph[SES_T1] != node_T0_graph[node_t] and check_edge[node_T0_graph[node_t]][node_T1_graph[SES_T1]] == 0)
            {
                add_edge(node_t, node_T1_graph[SES_T1], graph_of_shape);
                check_edge[node_T0_graph[node_t]][node_T1_graph[SES_T1]] = 1;
                parent_graph[node_t].push_back(node_T1_graph[SES_T1]);
            }

            if (node_T2_graph[SES_T2] != node_T0_graph[node_t] and check_edge[node_T0_graph[node_t]][node_T2_graph[SES_T2]] == 0)
            {
                add_edge(node_t, node_T2_graph[SES_T2], graph_of_shape);
                check_edge[node_T0_graph[node_t]][node_T2_graph[SES_T2]] = 1;
                parent_graph[node_t].push_back(node_T2_graph[SES_T2]);
            }
        }
    }




    // second tree
    for (int i = 0; i < representation_T1.size() ; i++)
    {
        if (!T[1].is_root(representation_T1[i]))
        {
            unsigned node_t = representation_T1[i];
            if (check_edge[node_T1_graph[node_t]][node_T1_graph[T[1].parent[node_t]]] ==0)
            {
                add_edge(node_T1_graph[node_t], node_T1_graph[T[1].parent[node_t]], graph_of_shape);
                check_edge[node_T1_graph[node_t]][node_T1_graph[T[1].parent[node_t]]] =1;
                parent_graph[node_T1_graph[node_t]].push_back(node_T1_graph[T[1].parent[node_t]]);

            }
            unsigned depth_temp = depth_graph[node_T1_graph[T[1].parent[node_t]]] + 1;
            if (depth_temp > depth_graph[node_T1_graph[node_t]])
                depth_graph[node_T1_graph[node_t]] = depth_temp;
            unsigned SES_T1 = SES[2][node_t];
            unsigned SES_T2 = SES[3][node_t];

            if (node_T0_graph[SES_T1] != node_T1_graph[node_t] and check_edge[node_T1_graph[node_t]][node_T0_graph[SES_T1]] == 0)
            {
                add_edge(node_T1_graph[node_t], node_T0_graph[SES_T1], graph_of_shape);
                check_edge[node_T1_graph[node_t]][node_T0_graph[SES_T1]] = 1;
                parent_graph[node_T1_graph[node_t]].push_back(node_T0_graph[SES_T1]);
            }

            if (node_T2_graph[SES_T2] != node_T1_graph[node_t] and check_edge[node_T1_graph[node_t]][node_T2_graph[SES_T2]] == 0)
            {
                add_edge(node_T1_graph[node_t], node_T2_graph[SES_T2], graph_of_shape);
                check_edge[node_T1_graph[node_t]][node_T2_graph[SES_T2]] = 1;
                parent_graph[node_T1_graph[node_t]].push_back(node_T2_graph[SES_T2]);
            }
        }

    }

    // third tree
    for (int i = 0; i < representation_T2.size() ; i++)
    {
        if (!T[2].is_root(representation_T2[i]))
        {
            unsigned node_t = representation_T2[i];
            if (check_edge[node_T2_graph[node_t]][node_T2_graph[T[2].parent[node_t]]] ==0)
            {
                add_edge(node_T2_graph[node_t], node_T2_graph[T[2].parent[node_t]], graph_of_shape);
                check_edge[node_T2_graph[node_t]][node_T2_graph[T[2].parent[node_t]]] =1;
                parent_graph[node_T2_graph[node_t]].push_back(node_T2_graph[T[2].parent[node_t]]);
            }
            unsigned depth_temp = depth_graph[node_T2_graph[T[2].parent[node_t]]] + 1;
            if (depth_temp > depth_graph[node_T2_graph[node_t]])
                depth_graph[node_T2_graph[node_t]] = depth_temp;
            unsigned SES_T1 = SES[4][node_t];
            unsigned SES_T2 = SES[5][node_t];

            if (node_T0_graph[SES_T1] != node_T2_graph[node_t] and check_edge[node_T2_graph[node_t]][node_T0_graph[SES_T1]] == 0)
            {
                add_edge(node_T2_graph[node_t], node_T0_graph[SES_T1], graph_of_shape);
                check_edge[node_T2_graph[node_t]][node_T0_graph[SES_T1]] = 1;
                parent_graph[node_T2_graph[node_t]].push_back(node_T0_graph[SES_T1]);
            }

            if (node_T1_graph[SES_T2] != node_T2_graph[node_t] and check_edge[node_T2_graph[node_t]][node_T1_graph[SES_T2]] == 0)
            {
                add_edge(node_T2_graph[node_t], node_T1_graph[SES_T2], graph_of_shape);
                check_edge[node_T2_graph[node_t]][node_T1_graph[SES_T2]] = 1;
                parent_graph[node_T2_graph[node_t]].push_back(node_T1_graph[SES_T2]);
            }
        }
    }


    //std::vector<unsigned> depth_node(max_graph+1);
    for(int i = 1; i < node_graph.size() ; i++)
    {
        //std::cout << node_graph[i] << std::endl;
        unsigned depth_temp = 0;
        for(int j = 0; j < parent_graph[node_graph[i]].size(); j++)
        {

            if (depth_temp < depth_graph[parent_graph[node_graph[i]][j]] + 1)
                depth_temp = depth_graph[parent_graph[node_graph[i]][j]] + 1;
        }
        depth_graph[node_graph[i]] = depth_temp;
        std::cout << node_graph[i]  << "   "  << depth_graph[node_graph[i]] << std::endl;
    }

//    for(int i = 0; i < node_graph.size() ; i++)
//    {
//        std::cout << node_graph[i] << "   "  << depth_graph[node_graph[i]] << std::endl;
//    }


//    // Reduction step
//    // remove unneccesary edges

    for(int i = 1; i < node_graph.size() ; i++)
    {
        std::cout << "node  "<<node_graph[i] << std::endl  ;
        unsigned max_depth = 0;
        unsigned index_maxdepth =0;
        for (int j = 0; j < parent_graph[node_graph[i]].size() ; j++)
        {
            if (depth_graph[parent_graph[node_graph[i]][j]] > max_depth)
            {
                max_depth = depth_graph[parent_graph[node_graph[i]][j]];
                index_maxdepth = j;
            }
        }

        for (int j = 0; j < parent_graph[node_graph[i]].size() ; j++)
        {
            if (j != index_maxdepth)
                remove_edge(node_graph[i], parent_graph[node_graph[i]][j], graph_of_shape);
        }
    }






////    // print edges
    std::cout << "edges(g) = ";
    IndexMap index_graph2 = get(vertex_index, graph_of_shape);
    graph_traits<UndirectedGraph>::edge_iterator ei1, ei_end1;
    for (boost::tie(ei1, ei_end1) = edges(graph_of_shape); ei1 != ei_end1; ++ei1)
        std::cout << "(" << index_graph2[source(*ei1, graph_of_shape)]
                  << "," << index_graph2[target(*ei1, graph_of_shape)] << ") ";
    std::cout << std::endl;







    ////////////////////////////////////////////////////////// Fill MAX TREE ///////////////////////////////////////
    //// build a tree from depth image


























//    typedef std::pair<uint8,uint8> pair_t;
//    std::vector<pair_t> minmax(finalNumberOfLabels1);
//    std::vector< std::vector<unsigned> > children(finalNumberOfLabels1);
//    std::vector<bool> enqueued(finalNumberOfLabels1);
//    std::fill(enqueued.begin(),enqueued.end(), false);
//    std::vector<uint8> dist1(finalNumberOfLabels1);

//    for (unsigned i = 0; i < T[2].S.size(); ++i)
//    {
//        unsigned sp = T[2].S[i];

//        minmax[sp] = pair_t(T[2].Ub[sp], T[2].Ub[sp]);
//        unsigned q = T[2].parent[sp];
//        children[q].push_back(sp);

//    }

//    unsigned start = label_border;



//    std::queue<unsigned> Qu;
//    Qu.push(start);
//    enqueued[start] = true;

//    while (not Qu.empty())
//    {
//        unsigned p = Qu.front();
//        Qu.pop();
//        dist1[p] = minmax[p].second - minmax[p].first;
//        //std::cout << p  << "   " <<int(dist1[p]) << std::endl;
//        // we do not process the nodes wrt increasing 'dist'...

//        unsigned par_p = T[2].parent[p];

//        if (not enqueued[par_p])
//        {
//           // std::cout << "par_p  "  <<par_p << std::endl;

//            if (par_p != p) // except for the root node
//              //mm(par_p).update_with(mm(p));
//            {
//                if (minmax[par_p].first > minmax[p].first )
//                    minmax[par_p].first = minmax[p].first;

//                if (minmax[par_p].second < minmax[p].second )
//                    minmax[par_p].second = minmax[p].second;
//            }

//            Qu.push(par_p);
//            enqueued[par_p] = true;
//        }

//        for (int t = 0; t < children[p].size() ; t++)
//        {
//            unsigned c = children[p][t];
//            if (not enqueued[c])
//            {
//                //std::cout << "c   " << c  << "  " <<int(minmax[c].first )   << "  "  << int (minmax[c].second)  << "   "<<int( T.Ub[p]) << std::endl;

//                  //mm(c).update_with(mm(p));
//                  if (minmax[c].first > minmax[p].first )
//                      minmax[c].first = minmax[p].first;

//                  if (minmax[c].second < minmax[p].second )
//                      minmax[c].second = minmax[p].second;

//                  Qu.push(c);
//                  enqueued[c] = true;
//            }

//        }
//    }



//    image2d<uint8_t>  dmap1(D);


//    for (int t = 0; t< finalNumberOfLabels; t++)
//    {
//        for(int e = 0; e < listpixel[t].size() ; e++)
//        {
//            unsigned py = listpixel[t][e]/width;
//            unsigned px = listpixel[t][e]%width;
//            dmap1(point2d(py,px)) = minmax[t].second - minmax[t].first;

//        }
//    }


//    io::imsave(dmap1, "vkl1.ppm");







//    double stop_s=clock();

//    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;


////    io::imsave(adjcMatrix, "vkl.ppm");
//    io::imsave(Seg_Img_gray, "dkm.ppm");
//    io::imsave(outlabels, "label.ppm");






}
