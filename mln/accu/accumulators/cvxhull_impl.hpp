#ifndef CONVEXHULL_HPP
# define CONVEXHULL_HPP

#include <mln/core/point.hpp>
#include <vector>
#include <cmath>

namespace mln
{

  /// \brief Compute angles from a reference point and retourne the one of maximum angle.
  /// \param pbegin Iterator to the begining of points array
  /// \param pend Iterator to the end of points array
  /// \param pend Iterator to the begining of output array to store angles
  /// \param ref reference point
  /// \return an iterator of the output array that points to its maximum
  template <typename InputIterator, typename OutputIterator>
  OutputIterator
  get_max_angles(InputIterator pbegin, InputIterator pend, point2d ref, OutputIterator out)
  {
    bool forward = (ref[0] <= (*pbegin)[0]);

    OutputIterator res = out;
    for (; pbegin != pend; ++pbegin, ++out)
      *out = std::atan( double((*pbegin)[1] - ref[1]) / ((*pbegin)[0] - ref[0]) );

    // FIXME: use soft (eps) float comparison to avoid multiple colinear segments

    if (forward) {
      for(OutputIterator it = res; it != out; ++it)
	if (*res <= *it) // non strict !
	  res = it;
    } else {
      for(OutputIterator it = res; it != out; ++it)
	if (*res < *it) // strict !
	  res = it;
    }

    return res;
  }



  ///
  /// \param points
  ///
  /// \pre \p points must be sorted along first dimension, then second dim.
  ///         (scan order)
  std::vector<point2d>
  convexhull(const std::vector<point2d>& points)
  {
    std::vector<point2d> cvx_hull;
    std::vector<float> angles(points.size());

    int i = 0;
    int n = points.size();
    point2d start, p;

    if (n == 0)
      return cvx_hull;

    start = points[0]; // upper left point
    cvx_hull.push_back(start);
    //std::cout << "UL: " << 0 << "," << start << std::endl;

    // Add colinear point
    i = 1;
    while (i < n && points[i][0] == start[0])
      ++i;
    if (i > 1 && i < n) {
      cvx_hull.push_back(points[i-1]);
      //std::cout << "insert: " << (i-1) << "," << points[i-1] << std::endl;
    }

    //std::cout << "UR: " << (i-1) << "," << points[i-1] << std::endl;

    // Go down
    {
      p = cvx_hull.back();
      while (i < n)
	{
	  auto it = get_max_angles(points.begin() + i, points.end(), p, angles.begin() + i);
	  unsigned pos = (it - angles.begin());
	  p = points[pos];
	  cvx_hull.push_back(p);
	  i = pos+1;
	}
    }

    //std::cout << "LR: " << (i-1) << "," << points[i-1] << std::endl;

    // Process colinear
    {
      while (i > 0 and points[i-1][0] == p[0])
	--i;
      if (i+1 != (n-1) and i >= 0)
	cvx_hull.push_back(points[i]);
      p = cvx_hull.back();
    }

    //std::cout << "LL: " << i << "," << points[i] << std::endl;
    // Go up
    // here i is past-the-end
    {
      while (i > 0)
	{
	  auto it = get_max_angles(points.begin(), points.begin() + i, p, angles.begin());
	  int pos = (it - angles.begin());
	  i = pos;
	  p = points[i];
	  cvx_hull.push_back(p);
	  //std::cout << "insert: " << i << "," << p << std::endl;
	}
    }

    //std::cout << "Loop: " << i << "," << points[i] << std::endl;

    return cvx_hull;
  }

}


#endif // ! CONVEXHULL_HPP
