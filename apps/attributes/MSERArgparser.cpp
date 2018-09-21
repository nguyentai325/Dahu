#include "MSERArgparser.hpp"
#include <boost/program_options.hpp>

MSERArgParser::MSERArgParser()
  : m_desc("MSER options")
{
  namespace po = boost::program_options;

  m_desc.add_options()
    ("mode,m", po::value<std::string>()->default_value("MSER_RATIO"),
     "Mode: (MSER_DIFF | MSER_RATIO | MSER_NORM) : \n"
     "MSER_DIFF: Area(parent) - Area(current)\n"
     "MSER_RATIO 0-1: Area(current)/Area(parent)\n"
     "MSER_NORM: Area(parent) / Area(current) - 1"
     )
    ("delta,d", po::value<float>()->default_value(20), "Number of gray levels in to accept a parent.")
    ("a0", po::value<float>()->default_value(5), "A0 Penalty term: P = a0.exp(a1 * area)")
    ("a1", po::value<float>()->default_value(0.0005), "A1 Penalty term: P = a0.exp(a1 * area)")
    ;
}

boost::program_options::options_description&
MSERArgParser::description()
{
  return m_desc;
}

