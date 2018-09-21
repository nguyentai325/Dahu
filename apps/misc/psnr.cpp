#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/accu/accumulators/sum.hpp>

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " input1 input2\n";
    std::exit(1);
  }


  image2d<rgb8> f, g;
  io::imread(argv[1], f);
  io::imread(argv[2], g);

  auto f_ = imcast<rgb<float>>(f);
  auto g_ = imcast<rgb<float>>(g);
  auto diff = imtransform(f_ - g_, [](rgb<float> x) -> double { return l2norm_sqr(x); });

  auto dims = f.domain().shape();
  double sum = accumulate(diff, accu::features::sum<double> ());
  double MSE = sum / (3 * dims[0] * dims[1]);
  std::cout << "MSE: " << MSE << "\n";

  if (MSE == 0)
    std::cout << "PSNR: " << -1 << "\n";
  else {
    float psnr = 20 * std::log10(255) - 10 * std::log10(MSE);
    std::cout << "PSNR: " << psnr << "\n";
  }
}
