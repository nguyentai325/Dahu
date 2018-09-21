#ifndef APPS_ATTRIBUTES_MSERARGPARSER_HPP
# define APPS_ATTRIBUTES_MSERARGPARSER_HPP

# include "argparser.hpp"
# include "MSER.hpp"
# include <boost/program_options/errors.hpp>
# include <boost/program_options/variables_map.hpp>
# include <apps/tos/topology.hpp>
# include <string>



class MSERArgParser : public AttributeArgParser
{
public:
  MSERArgParser();
  virtual boost::program_options::options_description& description();

  template <typename V, typename T>
  mln::image2d<float>
  run(const boost::program_options::variables_map& vm,
      const mln::image2d<V>& f,
      const mln::image2d<T>& K,
      const mln::image2d<unsigned>& parent,
      const std::vector<unsigned>& S);

private:
  boost::program_options::options_description m_desc;
};

/*****************************/
/** Implementation         ***/
/*****************************/

template <typename V, typename T>
mln::image2d<float>
MSERArgParser::run(const boost::program_options::variables_map& vm,
		   const mln::image2d<V>& f,
		   const mln::image2d<T>& K,
		   const mln::image2d<unsigned>& parent,
		   const std::vector<unsigned>& S)
{
  std::string smode = vm["mode"].as<std::string>();
  eMSER_attribute mode;
  if (smode == "MSER_DIFF")
    mode = MSER_DIFF;
  else if (smode == "MSER_RATIO")
    mode = MSER_RATIO;
  else if (smode == "MSER_NORM")
    mode = MSER_NORM;
  else
    throw boost::program_options::invalid_option_value("MSER mode");

  float a0 = vm["a0"].as<float>();
  float a1 = vm["a1"].as<float>();
  float delta = vm["delta"].as<float>();

  mln::image2d<float> amser = compute_MSER_attribute(f, K, parent, S, delta, mode);

  // Compute area attribute
  mln::image2d<unsigned> area;
  mln::resize(area, K).init(0);
  area[S[0]] = 1; // root handling
  for (int i = S.size()-1; i > 0; --i)
  {
    unsigned x = S[i];
    if (mln::K1::is_face_2(K.point_at_index(x)))
      area[x]++;
    area[parent[x]] += area[x];
  }

  auto attr = transform(imzip(amser, area),
			[a0, a1](const std::tuple<float,unsigned>& v) -> float{
			  return std::max(std::get<0>(v), a0 * std::exp(-a1 * std::get<1>(v)) );
			});
  return eval(attr);
}




#endif // ! APPS_ATTRIBUTES_MSERARGPARSER_HPP
