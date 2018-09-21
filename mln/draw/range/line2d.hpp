#ifndef LINE2D_HPP
# define LINE2D_HPP

# include <mln/core/iterator/iterator_base.hpp>
# include <mln/core/range/iterator_range.hpp>
# include <mln/core/point.hpp>

namespace mln
{

  namespace draw
  {

    struct line2d_iterator;
    typedef iterator_range<line2d_iterator> line2d_range;
    line2d_range line2d(const point2d& pbegin, const point2d& pend);

    /*******************/
    /** Implementation */
    /*******************/

    namespace internal
    {
      template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
      }
    }


    struct line2d_iterator : iterator_base<line2d_iterator, point2d, point2d>
    {
      line2d_iterator()
      {
      }

      line2d_iterator(const point2d& pbegin, const point2d& pend)
	: m_pbegin (pbegin), m_pend(pend)
      {
	point2d dp = pend - pbegin;
	m_octhan = false;
	if (std::abs(dp[1]) < std::abs(dp[0])) {
	  std::swap(dp[1], dp[0]);
	  m_octhan = true;
	}

	m_sy = internal::sgn(dp[0]);
	m_sx = internal::sgn(dp[1]);
	m_dy = std::abs(dp[0]);
	m_dx = std::abs(dp[1]);
	m_ddy = 2*m_dy;
	m_ddx = 2*m_dx;
	mln_assertion(m_dy <= m_dx);
      }


      void init()
      {
	if (m_octhan) {
	  m_p = {m_pbegin[1], m_pbegin[0]};
	} else {
	  m_p = m_pbegin;
	}
	m_e = m_ddy - m_dx;
	m_i = m_dx;
      }

      void next()
      {
	m_p[1] += m_sx;
	if (m_e > 0)
	  {
	    m_p[0] += m_sy;
	    m_e -= m_ddx;
	  }
	m_e += m_ddy;
	--m_i;
      }


      bool finished() const
      {
	return m_i < 0;
      }

      point2d dereference() const
      {
	if (m_octhan)
	  return point2d{m_p[1], m_p[0]};
	else
	  return m_p;
      }


    private:
      bool m_octhan;
      point2d m_pbegin;
      point2d m_pend;
      point2d m_p;
      int m_sx;
      int m_sy;
      int m_dx;
      int m_dy;
      int m_ddx;
      int m_ddy;
      int m_e;
      int m_i;
    };

    inline
    line2d_range
    line2d(const point2d& pbegin, const point2d& pend)
    {
      return line2d_range(line2d_iterator(pbegin, pend));
    }


  }

}
#endif // ! LINE2D_HPP
