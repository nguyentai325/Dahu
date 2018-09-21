#include <mln/core/image/image2d.hpp>
#include <mln/labeling/blobs.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/grays.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <mln/core/algorithm/transform.hpp>

int main(int argc, const char* argv[])
{
  using namespace mln;
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " input.pbm output.tiff" << std::endl;
    abort();
  }

  image2d<bool> ima;
  io::imread(argv[1], ima);

  trace::verbose = true;
  labeling::blobs(ima, c8, 0u);

  {
    image2d<unsigned> out;
    resize(out, ima).init(0u);
    labeling::impl::generic::blobs_boundcheck(ima, c8, 0u, out);

    io::imsave(transform(out, [](unsigned x) -> float { return x; }), argv[2]);
  }

}
