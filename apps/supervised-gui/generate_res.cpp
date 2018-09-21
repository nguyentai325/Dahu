#include <mln/core/image/image2d.hpp>
#include <mln/core/se/ball.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/morpho/structural/gradient.hpp>
#include <mln/core/algorithm/fill.hpp>
#include <mln/colors/literal.hpp>



int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " rec.ppm marked.ppm out.ppm\n"
      ;
    std::exit(1);
  }

  image2d<rgb8> F, M;

  io::imread(argv[1], F);
  io::imread(argv[2], M);


  se::ball2d b = se::make_ball2d(5.3);
  image2d<bool> M2;
  M2 = morpho::structural::external_gradient(M != rgb8(literal::zero), b,
                                             productorder_less<bool> (),
                                             [] (int x) -> bool { return x != 0; }
                                             );

  fill(F | M2, colors::literal::white);

  io::imsave(F, argv[3]);
}
