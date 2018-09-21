#ifndef APPS_ATTRIBUTES_ARGPARSER_HPP
# define APPS_ATTRIBUTES_ARGPARSER_HPP

# include <boost/program_options/options_description.hpp>


class AttributeArgParser
{
public:

  virtual boost::program_options::options_description& description() = 0;
};


#endif // ! APPS_ATTRIBUTES_ARGPARSER_HPP
