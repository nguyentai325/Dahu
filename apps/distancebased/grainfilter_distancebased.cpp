#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/colors/lab.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/algorithm/transform.hpp>

#include <mln/morpho/tos/tos.hpp>

#include <apps/tos/addborder.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/routines.hpp>
#include <apps/tos/topology.hpp>

#include <iostream>
#include <string>
#include <boost/format.hpp>

#include "compute_tree_distance.hpp"


void usage(char** argv)
{
  std::cout << argv[0] << " input output[wo extension] λ₁ [λ₂ [λ₃ ...] ]"
	    << std::endl;
  std::terminate();
}

int main(int argc, char** argv)
{
  if (argc < 4)
    usage(argv);


  const char* infname = argv[1];
  const char* outfname = argv[2];

  using namespace mln;

  image2d<rgb8> input;
  io::imread(infname, input);


  typedef lexicographicalorder_less<rgb8> Compare;
  Compare cmp;


  image2d<rgb8> Ori = addborder(input, cmp);
  image2d<rgb8> f = interpolate_k1(Ori);


  image2d<rgb8> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  auto dist = [](rgb8 x, rgb8 y) { return l1norm(x-y); };

  std::tie(K, parent, S) = morpho::ToSdistance(Ori, c4, point2d{0,0}, dist);

  // io::imprint(K);
  // io::imprint(parent);

  //io::imsave(K, "/tmp/K.tiff");

  auto A = attribute_compute(f, K, parent, S, accu::features::count<> (),
   			     K1::is_face_2);

  image2d<rgb8> tmp;
  image2d<rgb8> out = clone(f);
  resize(tmp, input);


  for (int i = 3; i < argc; ++i)
    {
      unsigned t = std::atoi(argv[i]);
      for (unsigned x: S)
	{
	  unsigned y = parent[x];
	  if (A[x] < t)
	    {
	      if (A[y] >= t)
		    out[x] = K[y];
	      else
		out[x] = out[y];
	    }
	}
      copy(out | sbox2d{out.domain().pmin + point2d{2,2}, out.domain().pmax - point2d{2,2}, {2,2} },
	   tmp);
      io::imsave(tmp, (boost::format("%s-%05i.tiff") % outfname % t).str().c_str());
    }
}
