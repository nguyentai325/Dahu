#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/paste.hpp>
#include <mln/core/algorithm/copy.hpp>
#include <mln/morpho/saturate.hpp>
#include <mln/labeling/blobs.hpp>
#include <mln/labeling/accumulate.hpp>
#include <mln/accu/accumulators/count.hpp>
#include <mln/colors/literal.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char** argv)
{
  using namespace mln;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("no-sat", "do not saturate")
    ;

  po::options_description hidden;
  desc.add_options()
    ("mask.pbm", po::value<std::string> ())
    ("input.ppm", po::value<std::string> ())
    ("out.ppm", po::value<std::string> ())
    ;

  po::positional_options_description pd;
  pd.add("mask.pbm", 1)
    .add("input.ppm", 1)
    .add("out.ppm", 1);

  po::options_description allopt;
  allopt.add(desc).add(hidden);

  po::variables_map vm;

  try {
    po::store(po::basic_command_line_parser<char>(argc, argv).options(allopt).positional(pd).run(), vm);
    po::notify(vm);
  }
  catch (...) {
    std::cout << "Usage: " << argv[0] << " [--no-sat] mask.pbm input.ppm output.ppm\n"
              << "Keep only the largest CC and saturate\n";
    std::exit(1);
  }


  image2d<bool> mask;
  image2d<rgb8> f;

  io::imread(vm["mask.pbm"].as<std::string> (), mask);
  io::imread(vm["input.ppm"].as<std::string> (), f);
  uint16 nlabel;
  image2d<uint16> lbl;

  mask = eval(lnot(mask));
  std::tie(lbl, nlabel) = labeling::blobs(mask, c8, (uint16)0);

  auto res = labeling::p_accumulate(lbl, nlabel, accu::features::count<> ());

  res[0] = 0;
  std::vector<int> idx(nlabel);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&res](int i, int j) { return res[i] > res[j]; });

  int obj = idx[0];
  fill(mask, false);
  for (int i = 0; (i < nlabel) and (res[idx[i]] >= 0.5 * res[obj]); ++i)
    {
      std::cout << "Objet: " << idx[i] << " / " << nlabel << std::endl;
      fill((mask | (lbl == idx[i])), true);
    }

  if (vm.count("no-sat") == 0)
    {
      box2d dom2 = mask.domain();
      dom2.pmin -= 1;
      dom2.pmax += 1;
      image2d<bool> mask2(dom2, 3, false);
      paste(mask, mask2);
      mask2 = morpho::saturate(mask2, c8, point2d{-1,-1});
      copy(mask2 | mask.domain(), mask);
    }

  fill(f | lnot(mask), colors::literal::black);

       io::imsave(f, vm["out.ppm"].as<std::string> ());
}
