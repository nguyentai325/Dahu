#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/grays.hpp>

#include <mln/core/algorithm/transform.hpp>

#include <mln/morpho/tos/tos.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <apps/tos/addborder.hpp>

#include <apps/attributes/MSERArgparser.hpp>
#include <apps/attributes/meaningfullnessArgparser.hpp>

#include <boost/program_options.hpp>
#include "extinction.hpp"
#include "saliency.hpp"

void usage(char** argv)
{
  std::cout << argv[0] << " method input.tiff output.tiff" << std::endl
	    << "method = (mser|meaningfullness)" << std::endl;
  std::terminate();
}

int main(int argc, char** argv)
{
  using namespace mln;
  namespace po = boost::program_options;


  MSERArgParser mseroptions;
  MeaningFullNessArgParser meaningfullness_options;
  po::options_description mserdesc = mseroptions.description();
  po::options_description meaningfullness_desc = meaningfullness_options.description();

  po::options_description hidden;
  hidden.add_options()
    ("method", po::value<std::string>()->required(), "Method (mser|meaningfullness)")
    ("input", po::value<std::string>()->required(), "Input")
    ("output", po::value<std::string>()->required(), "Output (tiff)")
    ("ignore", po::value<std::vector<std::string> >(), "")
    ;

  po::positional_options_description pd;
  pd.add("method", 1);
  pd.add("input",  1);
  pd.add("output", 1);
  pd.add("ignore", -1);

  po::options_description shown("Allowed options");
  shown.add(mserdesc).add(meaningfullness_desc);

  po::options_description all;
  all.add(mserdesc).add(meaningfullness_desc).add(hidden);


  try {
    po::variables_map vm;
    po::command_line_parser parser(argc, argv);
    po::store(parser.options(hidden).positional(pd).allow_unregistered().run(), vm);
    vm.notify();

    // Load the image and compute the tos
    std::string ifname = vm["input"].as<std::string>();
    std::string ofname = vm["output"].as<std::string>();

    image2d<uint8> ima_, ima;
    io::imread(ifname.c_str(), ima_);
    ima = addborder(ima_);

    typedef UInt<9> V;
    typedef image2d<V> I;

    I f = transform(ima, [](uint8 v) -> V { return v * 2; });

    image2d<V> K;
    image2d<unsigned> parent;
    std::vector<unsigned> S;
    std::tie(K, parent, S) = morpho::ToS(f, c4);


    // Run the method
    std::string method = vm["method"].as<std::string>();
    image2d<float> attr;
    if (method == "mser")
      {
	po::variables_map vm;
	po::command_line_parser parser(argc, argv);
	po::store(parser.options(hidden.add(mserdesc)).positional(pd).run(), vm);
	vm.notify();
	attr = mseroptions.run(vm, K, K, parent, S);
      }
    else if (method == "meaningfullness")
      {
	po::variables_map vm;
	po::command_line_parser parser(argc, argv);
	po::store(parser.options(hidden.add(meaningfullness_desc)).positional(pd).run(), vm);
	vm.notify();
	attr = meaningfullness_options.run(vm, ima, K, parent, S);
      }
    else
      {
	throw po::invalid_option_value("MSER mode");
      }

    image2d<float> extmap = extinction(attr, K, parent, S, std::less<float> ());
    //collapse_zero_nodes(extmap, K, parent, S);
    image2d<float> s = saliencymap(extmap, K, parent, S);

    {
      image2d<float> tmp = clone(extmap);
      for (unsigned x: S)
	{
	  if (K[x] == K[parent[x]] or tmp[x] == 0)
	    tmp[x] = tmp[parent[x]];
	}
      io::imsave(tmp, "tmp.tiff");
    }


    io::imsave(s, ofname.c_str());

  } catch (po::error e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: " << argv[0] << " method input output" << std::endl;
    std::cerr << shown;
  }
}
