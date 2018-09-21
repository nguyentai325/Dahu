#include "libtos.hpp"

template
void
grainfilter_inplace<mln::uint8>(mln::image2d<mln::uint8>& K,
				mln::image2d<unsigned>& parent,
				const std::vector<unsigned>& S,
				float threshold);

template
void
grainfilter_inplace< mln::UInt<9> >(mln::image2d< mln::UInt<9> >& K,
				    mln::image2d<unsigned>& parent,
				    const std::vector<unsigned>& S,
				    float threshold);

template
void
grainfilter_inplace< unsigned >(mln::image2d< unsigned >& K,
				mln::image2d<unsigned>& parent,
				const std::vector<unsigned>& S,
				float threshold);
namespace mln
{

  namespace morpho
  {

    // template
    // std::tuple< image2d< UInt<9> >, image2d<unsigned>, std::vector<unsigned> >
    // ToS< image2d< UInt<9> >, c4_t, std::less< UInt<9> > >
    //   (const Image< mln::image2d<UInt<9> > >& ima,
    //    const c4_t& nbh,
    //    const std::less< UInt<9> >& cmp);
  }

}
