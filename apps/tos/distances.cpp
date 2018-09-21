# include <mln/core/image/image2d.hpp>
# include <mln/io/imread.hpp>

# include <mln/core/algorithm/accumulate.hpp>
# include <mln/accu/accumulators/max.hpp>

# include <apps/tos/colorToSGrad.hpp>
# include <apps/tos/addborder.hpp>
# include <apps/tos/utils.hpp>
# include <apps/tos/shape_similarity.hpp>

#include <mln/io/imsave.hpp>
#include <map>

void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input(color) segmentatio.tiff" << std::endl
	    << std::endl;
  std::terminate();
}

namespace mln
{

  template <typename V>
  image2d<V>
  add_border_and_interp_nn_x2(const image2d<V>& ima, const V& border_val)
  {
    box2d dom = ima.domain();
    dom.pmax = (dom.pmax + 2) * 2 - 1;

    image2d<V> out(dom, 3, border_val);
    mln_foreach(const point2d& p, ima.domain()) {
      point2d q = 2 * (p + point2d{1,1});
      out(q) = ima(p);
      out.at(q + point2d{0,1}) = ima(p);
      out.at(q + point2d{1,0}) = ima(p);
      out.at(q + point2d{1,1}) = ima(p);
    }
    return out;
  }


  template <typename V>
  image2d<uint8>
  labelize_uniq(const image2d<V>& ima)
  {
    std::map<V, uint8> map;
    int nlabel= 0;

    image2d<uint8> out;
    resize(out, ima);
    mln_pixter(pin, pout, ima, out);
    mln_forall(pin, pout) {
      auto res = map.insert( std::make_pair(pin->val(), nlabel) );
      if (res.second)
	nlabel++;
      pout->val() = res.first->second;
    }

    return out;
  }

}


int main(int argc, char** argv)
{
  if (argc < 3)
    usage(argv);

  using namespace mln;
  image2d<rgb8> ima;
  image2d<uint8> ref;
  io::imread(argv[1], ima);
  io::imread(argv[2], ref);

  image2d<rgb8> ima_ = addborder(ima, lexicographicalorder_less<rgb8>() );
  image2d<uint8> f = add_border_and_interp_nn_x2(ref, value_traits<uint8>::max());

  image2d<unsigned> K, K_;
  std::vector<unsigned> S, S_;
  image2d<unsigned> parent, parent_;
  colorToSGrad(ima, K_, parent_, S_);

  std::tie(K, parent, S) = remove_01face_from_tos(K_, parent_, S_);

  uint8 nlabel = accumulate(ref, accu::features::max<>() );
  coveryRate(f, nlabel, K, parent, S);


  // io::imsave(labelize_uniq(f), "labels.tiff");
  // io::imsave(labelize_uniq(K), "K.tiff");
}
