#include <mln/io/imread.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/morpho/component_tree/accumulate.hpp>
#include <mln/accu/accumulators/accu_if.hpp>
#include <mln/accu/accumulators/count.hpp>
#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/topology.hpp>
#include "attributes.hpp"
#include "cMeaningFullNess.hpp"



int main(int argc, char** argv)
{
  using namespace mln;

  const char* usage =
    "Use the energy E(Γ) = α V₁(Γ) + βV₂(Γ) + V₃(Γ)\n"
    "\twith\n"
    "\t     ε: The radius of the ball when considering the contextual energy (3-10).\n"
    "\t     V₁(Γ) = (ExternalVar(Γ,ε) + InternalVar(Γ,ε)) / Var(Γ,ε)\n"
    "\t     V₂(Γ) = Mean Curvature(∂Γ)\n"
    "\t     V₂(Γ) = exp(-γ |∂Γ|/|∂Ω|)\n"
    "\t     Consider α=1, β=1, γ=10-100\n"
    "\t     (you can set one of these parameters to 0 to unactive the term).";

  po::options_description desc("MSER Options");
  desc.add_options()
    ("params,p", po::value< std::vector<float> >()->multitoken()->required(), "ε α β γ")
    ;

  po::variables_map vm = process_cmdline(argc, argv, desc, usage);
  std::vector<float> params = vm["params"].as< std::vector<float> >();

  if (params.size() != 4) {
    std::cerr << usage;
    std::exit(1);
  }

  int eps = params[0];
  float alpha = params[1], beta = params[2], gamma = params[3];
  tree_t tree = preprocess(vm);

  typedef rgb8 V;

  image2d<V> f;
  io::imread(vm["input_path"].as<std::string>(), f, true);

  image2d<V> F = Kadjust_to(f, tree._get_data()->m_pmap.domain());
  f = unimmerse_k1(F);

  // Curvature
  image2d<float> C = compute_curvature(transform(f, [](V x) { return l1norm(x); }));


  // V1
  typedef std::conditional< std::is_scalar<V>::value, double,
                            vec<double, value_traits<V>::ndim> >::type SumType;
  auto amap = compute_regional_energy(tree, F,
                                      accu::features::variance<SumType> (),
                                      eps);
  amap[tree.get_root()] = 0;

  // V2
  auto amap2 = compute_attribute_on_contour(tree, C,
                                            accu::features::count<unsigned> () &
                                            accu::features::mean<double> ());


  // Compute energy
  property_map<tree_t, float> energy(tree);
  {
    float domlength = tree._get_data()->m_pmap.domain().size();
    mln_foreach(auto x, tree.nodes())
      {
        float v1 = amap[x][2] > 0 ? (alpha * (amap[x][0]+amap[x][1])/amap[x][2]) : 0;
        float v2 = beta * accu::extractor::mean(amap2[x]);
        float v3 = std::exp(-gamma * accu::extractor::count(amap2[x]) / domlength);
        energy[x] = v1+v2+v3;
      }
  }
  energy[tree.get_root()] = 0;


  auto smap = postprocess(vm, tree, energy);


  const char* names[] = {"VIN", "VOUT", "VTOTAL", "V1", "Cont. Length",
                         "Curv. Avg", "Energy", "Extinction" };

  typedef vec<double, 3> vec3d;
  auto vin = make_functional_property_map<tree_t::node_type>([&amap](tree_t::node_type x) { return amap[x][0]; });
  auto vout = make_functional_property_map<tree_t::node_type>([&amap](tree_t::node_type x) { return amap[x][1]; });
  auto vtotal = make_functional_property_map<tree_t::node_type>([&amap](tree_t::node_type x) { return amap[x][2]; });
  auto v1 = make_functional_property_map<tree_t::node_type>([&amap](tree_t::node_type x) {
      return amap[x][2] == 0 ? 0 : (amap[x][0]+amap[x][1]) / amap[x][2];
    });

  auto cont_length = make_functional_property_map<tree_t::node_type>([&amap2](tree_t::node_type x) {
      return accu::extractor::count(amap2[x]); });
  auto cont_curv = make_functional_property_map<tree_t::node_type>([&amap2](tree_t::node_type x) {
      return accu::extractor::mean(amap2[x]); });

  std::function<float(tree_t::node_type)> attrs[] = {
    _as_fun(vin),
    _as_fun(vout),
    _as_fun(vtotal),
    _as_fun(v1),
    _as_fun(cont_length),
    _as_fun(cont_curv),
    _as_fun(energy),
    _as_fun(smap)
  };

  export_(vm, tree, smap, names, attrs, sizeof(names) / sizeof(char*));
}
