#ifndef MLN_CORE_SE_BALL_HPP
# define MLN_CORE_SE_BALL_HPP

# include <mln/core/neighborhood/dyn_neighborhood.hpp>
# include <vector>

namespace mln
{
  namespace se
  {

    /// Ball structuring element.
    ///
    /// A ball of radius \$r\$ stands is centered in (0,0) and has all the
    /// points whose (euclidean) distance is lower than \$r\$
    struct ball2d;

    /// \brief Helper function to make a Ball
    ball2d make_ball2d(float r);

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    struct ball2d : dyn_neighborhood_base<
      std::vector<point2d>,
      dynamic_neighborhood_tag,
      ball2d>
    {
      using is_incremental = std::true_type;

      using dec_type = dyn_neighborhood<std::vector<point2d>,
                                        dynamic_neighborhood_tag>;
      using inc_type = dyn_neighborhood<std::vector<point2d>,
                                        dynamic_neighborhood_tag>;

      ball2d() = default;

    private:
      ball2d(const std::vector<point2d>& dpts,
             const std::vector<point2d>& vdec,
             const std::vector<point2d>& vinc)
        : dyn_neighborhood_base(dpts),
          m_dec(vdec),
          m_inc(vinc)
      {
      }


    public:
      static
      ball2d
      make(float r_)
      {
        typedef point2d::value_type P;
        mln_precondition(r_ >= 0);

        int k = r_;
        int d = 2*k+1;
        float r2 = r_ * r_;

        std::vector<point2d> vinc, vdec, dpoints;
        dpoints.reserve(d*d);
        vinc.reserve(d);
        vdec.reserve(d);

        for (int i = -k; i <= k; ++i)
          for (int j = -k; j <= k; ++j)
            if (sqr(i) + sqr(j) <= r2) {
              point2d p = {(P)i, (P)j};
              dpoints.push_back(p);
            }

        {
          int n = dpoints.size();
          for (int i = 0; i < n; i++)
            {
              point2d p = dpoints[i];
              vdec.push_back(p + point2d{0,-1}); // before begin of the line
              while (i+1 < n and dpoints[i][0] == dpoints[i+1][0])
                ++i;

              vinc.push_back(dpoints[i]); // last element of the line
            }
        }
        return ball2d(dpoints, vdec, vinc);
      }

      inc_type
      inc() const
      {
        return m_inc;
      }

      dec_type
      dec() const
      {
        return m_dec;
      }

    private:
      dec_type m_dec;
      inc_type m_inc;
    };

    inline
    ball2d
    make_ball2d(float r)
    {
      return ball2d::make(r);
    }

  } // end of namespace mln::se

} // end of namespace mln

#endif //!MLN_CORE_SE_BALL_HPP
