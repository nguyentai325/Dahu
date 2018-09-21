#include <mln/core/image/image2d.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/morpho/component_tree/io.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/morpho/component_tree/filtering.hpp>
#include <mln/morpho/component_tree/reconstruction.hpp>
#include <mln/accu/accumulator.hpp>
#include <mln/accu/accumulators/accu_if.hpp>
#include <mln/accu/accumulators/min.hpp>
#include <mln/accu/accumulators/max.hpp>
#include <mln/accu/accumulators/moment_of_inertia.hpp>
#include <mln/colors/literal.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include <apps/tos/croutines.hpp>
#include <iostream>
#include <boost/format.hpp>
typedef mln::morpho::component_tree<unsigned, mln::image2d<unsigned> > tree_t;

namespace mln
{


  template <class V, class VMap>
  std::vector< image2d<float> >
  profile(const tree_t& tree, const image2d<V>& F, const VMap& vmap, std::vector<float>& lambdas,
          std::vector< image2d<typename VMap::value_type> >& res)
  {
    typedef accu::accumulators::moment_of_inertia<point2df, vec2f> myaccu_t;

    auto amap = morpho::paccumulate(tree, F, myaccu_t ());
    typedef typename VMap::value_type W;

    float amap_max;
    {
      auto tmp = make_image(tree, amap);
      amap_max = accumulate(tmp, accu::features::max<>());

      for (float x = 0.1; x < amap_max; x += 0.1)
        lambdas.push_back(x);
    }

    res.reserve(lambdas.size()+1);
    {
      image2d<float> out;
      resize(out, F);
      morpho::reconstruction(tree, amap, out);
      io::imsave(out, "amap.tiff");
    }

    {
      image2d<W> out;
      resize(out, F);

      morpho::reconstruction(tree, vmap, out);
      //io::imsave(imcast<rgb16>(out), "rec_0.00.tiff");
      res.push_back(std::move(out));

      for (float t : lambdas)
        {
          image2d<W> out;
          resize(out, F);

          auto crit = make_functional_property_map<tree_t::node_type>([amap, t](tree_t::node_type x) { return amap[x] >= t; });
          morpho::filter_direct_and_reconstruct(tree, crit, vmap, out);

          //io::imsave(imcast<rgb16>(out), (boost::format("rec_%.02f.tiff") % t).str());
          res.push_back(std::move(out));
        }
    }

    std::vector< image2d<float> > out;
    out.resize(lambdas.size());

    for (unsigned i = 0; i < lambdas.size(); ++i) {
      out[i] = eval(imcast<float>(imtransform(res[i+1] - res[i], functional::l2norm_t<> ())));
    }

    return out;
  }

  template <typename T>
  struct mymean : Accumulator< mymean<T> >
  {
    typedef accu::accumulators::mean<T> accu_t;
    typedef void                         argument_type;
    typedef typename accu_t::result_type result_type;

    mymean() = default;

    void init() {
      m_accu.init();
    }

    template <typename Pix>
    void take(const Pix& pix) {
      if (K1::is_face_2(pix.point())) {
        m_accu.take(pix.val());
      }
    }

    void take(const mymean& other) {
      m_accu.take(other.m_accu);
    }


    result_type to_result() const {
      return m_accu.to_result();
    }

  private:
    accu_t m_accu;
  };

}
