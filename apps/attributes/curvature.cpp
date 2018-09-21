#include "curvature.hpp"

namespace mln
{

  template
  image2d<float> compute_curvature<uint8>(const image2d<uint8>&);

  template
  image2d<float> compute_curvature<uint16>(const image2d<uint16>&);
}
