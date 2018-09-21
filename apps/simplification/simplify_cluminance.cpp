#include <iostream>


#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/colors.hpp>

#include <mln/core/algorithm/transform.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <apps/tos/set_mean_on_nodes.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>



#include <apps/tos/addborder.hpp>
#include "simplify.hpp"

void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input(gray) eps output.tiff [tolerance]" << std::endl
	    << "Perform a simplification of the ToS by removing overlapping"
	    << "level lines in radius of lambda pixels" << std::endl
	    << "Param:" << std::endl
	    << "eps: size of the dilation" << std::endl
	    << "tolerance: Allow shapes that provides at least n% of non-overlap" << std::endl
	    << "           (0 <= n <= 1): 0 means all shapes, 1 means strictly included shapes." << std::endl
    //	    << "grainsize: minimal grain size (default: 0)" << std::endl
    //	    << "areafactor: avoid sur-simplification by limiting node desactivation (default: 0.0)" << std::endl
    //	    << "\t a shape S1 can only be desactived by a shape S2 if area(S1) * areafactor < area(S2)"
	    << std::endl
    ;
  std::terminate();
}


int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 4)
    usage(argv);

  int eps = std::atoi(argv[2]);
  float tol = argc >= 5 ? std::atof(argv[4]) : 1.0;


  image2d<rgb8> ori, ori_;
  io::imread(argv[1], ori_);
  ori = addborder(ori_, lexicographicalorder_less<rgb8> ());

  image2d<uint8> ima_, ima;
  ima = transform(ori, [](const rgb8& v) -> uint8 { return (v[0] + v[1] + v[2]) / 3; });


  typedef UInt<9> V;
  typedef image2d<V> I;

  I f = transform(ima, [](uint8 v) -> V { return v * 2; });

  image2d<V> K;
  image2d<unsigned> parent;
  std::vector<unsigned> S;

  std::tie(K, parent, S) = morpho::ToS(f, c4);


  //image2d<uint8> simp = simplify_bottom_up(ima, K, parent, S, lambda, grainsize, areafactor);
  image2d<rgb8>  mean = set_mean_on_node2(immerse_k1(ori), K, S, parent, K1::is_face_2);
  image2d<rgb8>  simp2 = simplify_top_down_tolerance(unimmerse_k1(mean), K, parent, S, eps, tol);

  //io::imsave(simp, argv[3]);
  io::imsave(simp2, argv[3]);
}
