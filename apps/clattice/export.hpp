#ifndef EXPORT_HPP
# define EXPORT_HPP

# include <boost/format.hpp>
# include <boost/numeric/ublas/matrix.hpp>

# include <mln/io/png/imsave.hpp>
# include <fstream>
# include "shape.hpp"



enum e_output_mode {
  FMT_CLUSTER = 0x01, // Outout shape as a cluster where each point is a node
  FMT_IMAGE = 0x02,   // Output shape as an image
  FMT_SUMMARY = 0x04  // Only summary (i.e for cluster only representant)
};

/// \brief Export a shape in dot as a cluster containing all points
/// of the shape as node.
template <typename shape_t>
std::string
export_shape_as_cluster(int id, const shape_t& shp, std::ostream& dotsream,
			mln::image2d<bool>& deja_vu,
			bool summary = true);


/// \brief Export a shape in dot as node containing the image of the
/// shape.
///
/// Note: the image is exported as pbm (shape_id.pbm) and has to be
/// converted to png before dottififcation
template <typename shape_t>
std::string
export_shape_as_image(int id, const shape_t& shp,
		      const mln::box2d& domain,
		      std::ostream& dotsream,
		      const std::string& dirname = "./");


template <typename shape_t, typename V>
void
save_graph(const boost::numeric::ublas::matrix<bool>& mat,
	   const std::vector<shape_t>& shapes,
	   const mln::image2d<V>& ima,
	   const std::string& filename);


/**************************/
/** Implementation        */
/**************************/


template <typename shape_t>
std::string
export_shape_as_cluster(int id,
			const shape_t& shp,
			std::ostream& dotsream,
			mln::image2d<bool>& deja_vu,
			bool summary)
{
  using namespace mln;

  // New cluster
  dotsream << "  subgraph cluster" << id << "{" << std::endl
	   << "label=\"" << shp << "\";" << std::endl;

  // Add shape points to the cluster
  bool has_repr = false;
  mln::point2d repr;
  for (const mln::point2d& p : shp.pset())
    if (!deja_vu(p)) {
      if (!summary or !has_repr)
	dotsream << '"' << p << '"' << " [style=filled, color=\".3 .7 1.\"]" << std::endl;
      if (!has_repr) {
	repr = p;
	has_repr = true;
      }
    }

  // This node does not provide any new points
  // We add a fake point for visualization
  std::string str_repr;
  if (!has_repr) {
      str_repr = (boost::format("empty_%i") % id).str();
      dotsream << str_repr << " [shape=point];" << std::endl;
  } else {
    str_repr = (boost::format("%s") % repr).str();
  }

  // Link to representant
  bool is_repr = true;
  for (const mln::point2d& p : shp.pset())
    {
      if (!deja_vu(p) and !summary)
	{
	  if (is_repr)
	    is_repr = false;
	  else
	    dotsream << '"' << p << "\" -> "
		     << '"' << str_repr << "\";" << std::endl;
	}
      deja_vu(p) = true;
    }

  dotsream << "  }" << std::endl;
  return str_repr;
}

template <typename shape_t>
std::string
export_shape_as_image(int id,
		      const shape_t& shp,
		      const mln::box2d& domain,
		      std::ostream& dotsream,
		      const std::string& dirname
		      )
{
  // New cluster
  dotsream << "  subgraph cluster" << id << "{" << std::endl
	   << "label=\"" << shp << "\";" << std::endl;

  dotsream << "   shape_" << id
	   << "[shape=rectangle,label=\"\","
	   << "image=\"" << dirname << "/" << "shape_" << id << ".png\", imagescale=true];"
	   << std::endl
	   << "  }" << std::endl;


  std::string repr = (boost::format("shape_%i") % id).str();

  mln::image2d<bool> bin(domain, 0, false);
  for (const mln::point2d& p: shp.pset())
    bin(p) = true;
  mln::io::png::imsave(bin, (dirname + "/" + repr + ".png").c_str());
  return repr;
}


template <typename shape_t, typename V>
void
save_graph(const boost::numeric::ublas::matrix<bool>& mat,
	   const std::vector<shape_t>& shapes,
	   const mln::image2d<V>& ima,
	   const std::string& filename,
	   const std::string& dirname)
{
  std::ofstream file(filename);
  mln::image2d<bool> deja_vu;
  mln::resize(deja_vu, ima).init(false);


  file << "digraph G {" << std::endl
       << "  compound = true;" << std::endl;


  std::vector<std::string> repr(shapes.size());

  int mode = FMT_IMAGE;
  for (unsigned i = 0; i < shapes.size(); ++i)
    {
      if (mode & FMT_CLUSTER)
	repr[i] = export_shape_as_cluster(i, shapes[i], file, deja_vu, mode & FMT_SUMMARY);
      else if (mode & FMT_IMAGE)
	repr[i] = export_shape_as_image(i, shapes[i], ima.domain(), file, dirname);
    }

  // Link groups
  int sz = shapes.size();
  for (int i = 0; i < sz; ++i)
    for (int j = i+1; j < sz; ++j)
      if (mat(i,j))
	file << "  \"" << repr[i] << "\" -> \"" << repr[j] << "\";" << std::endl;
  file << "}";
}



#endif // ! EXPORT_HPP
