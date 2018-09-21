#include "meaningfullnessArgparser.hpp"
#include <boost/program_options.hpp>

MeaningFullNessArgParser::MeaningFullNessArgParser()
  : m_desc("MeaningFullNess Attribute: α.Ei + Ee + Ec")
{
  namespace po = boost::program_options;

  m_desc.add_options()
    ("alpha", po::value<float>()->default_value(2), "Internal energy coef.")
    ("a0", po::value<float>()->default_value(1), "A0 Penalty term: Ec = a0.exp(a1 * area)")
    ("a1", po::value<float>()->default_value(0.005), "A1 Penalty term: Ec = a0.exp(a1 * area)")
    ("eps", po::value<int>()->default_value(5), "Radius for external/internal region")
    ;
}

boost::program_options::options_description&
MeaningFullNessArgParser::description()
{
  return m_desc;
}
