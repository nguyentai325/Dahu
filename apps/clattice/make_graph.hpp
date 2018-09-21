#ifndef MAKE_GRAPH_HPP
# define MAKE_GRAPH_HPP

 #define BOOST_UBLAS_TYPE_CHECK 0

# include <boost/numeric/ublas/matrix.hpp>
# include <boost/numeric/ublas/triangular.hpp>
# include <boost/numeric/ublas/matrix_proxy.hpp>
# include <boost/numeric/ublas/io.hpp>
# include <utility>

namespace mln
{

  // \param shapes is the set of (non redoundant) shapes
  // \param pt_cmp is a function to compare point of immersion space
  //
  // This algorithm proceeds as follow:
  // 1. It first sort shapes
  // 2. It compute the transitive reduction of the inclusion relation
  //    and return its Hass diagram as a an adjancy matrix
  // \returns (C,R) closure matrix / reduction matrix
  template <typename shape_t>
  std::pair<
    boost::numeric::ublas::triangular_matrix<bool, boost::numeric::ublas::upper>,
    boost::numeric::ublas::triangular_matrix<bool, boost::numeric::ublas::upper> >
  graph_transitive_reduction(std::vector<shape_t>& vs);



  /******************/
  /* Implementation */
  /******************/


  // The algorithm process
  template <typename shape_t>
  std::pair<
    boost::numeric::ublas::triangular_matrix<bool, boost::numeric::ublas::upper>,
    boost::numeric::ublas::triangular_matrix<bool, boost::numeric::ublas::upper> >
  graph_transitive_reduction  (std::vector<shape_t>& vs)
  {

    // mln_precondition( pt_cmp(point2d(0,0), point2d(0,1)) );
    // mln_precondition( pt_cmp(point2d(0,0), point2d(1,0)) );
    // mln_precondition( not pt_cmp(point2d(0,1), point2d(0,0)) );
    // mln_precondition( not pt_cmp(point2d(1,0), point2d(0,0)) );


    // By sorting, we get the property:
    //   A C B => Ord(A) < Ord(B)
    int s = vs.size();
    std::sort(vs.begin(), vs.end(), [](const shape_t& lhs, const shape_t& rhs)
	      { return lhs.size() < rhs.size(); });

    using namespace boost::numeric;
    ublas::triangular_matrix<bool, ublas::upper> R(s,s); // matrice of reduction
    ublas::triangular_matrix<bool, ublas::upper> C(s,s); // matrice of closure
    std::fill(C.data().begin(), C.data().end(), false);
    std::fill(R.data().begin(), R.data().end(), false);

    // 2. Make the matrice of inclusions
    // We are going to compute both transitive closure and transitive reduction
    for (int i = s-1; i >= 0; --i)
      for (int j = i+1; j < s; ++j)
	if (not C(i,j))
	  {
	    //mln_invariant(vs[i] != vs[j]); // Invariant very costly !
	    bool inc = is_shape_included(vs[i], vs[j]);
	    if (inc)
	      {
		R(i,j) = C(i,j) = true;
		ublas::subrange(C, i,i+1, j+1,s) += ublas::subrange(C, j,j+1, j+1,s); // Closure O(n)
	      }
	  }

    //std::cout << "Closure:" << std::endl << C << std::endl;
    //std::cout << "Reduction:" << std::endl << R << std::endl;

    return std::make_pair(C,R);
  }

}


#endif // ! MAKE_GRAPH_HPP
