#ifndef MLN_CONTRIB_MEANSHIFT_MEANSHIFT_HPP
# define MLN_CONTRIB_MEANSHIFT_MEANSHIFT_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/win2d.hpp>

namespace mln
{
  namespace contrib
  {


    template <class V>
    image2d<V>
    meanshift(const image2d<V>& f,
              float hs,
              float hr)
    {
      int SR = 5; // Spatial window radius
      int CR = 10; // Color window radius
      int NITER = 30; // Maximal number of iteration
      float eps = 0.1;  //
      double hs2 = hs*hs;
      double hr2 = hr*hr;
      double eps2 = eps * eps;

      typedef vec<double, value_traits<V>::ndim> value_t;
      typedef vec<double, 2> site_t;

      image2d<V> out;
      resize(out, f);

      rect2d win = make_rectangle2d(2*SR+1, 2*SR+1);

      point2d x;
      mln_iter(p__, f.domain());
      mln_iter(q__, win(x));

      auto g = [](double x) -> double { return std::exp(-x);  };
      mln_foreach(point2d p, p__)
        {
          site_t py = p.as_vec();
          value_t vy = f(p).as_vec();

          bool stop = false;
          for (int i = 0; i <= NITER and (not stop); ++i)
            {
              site_t s1 = literal::zero;
              value_t s2 = literal::zero;
              double s3 = 0;
              x = (point2d)py;
              mln_foreach(point2d q, q__)
                {
                  if (f.domain().has(q))
                    {
                      double d0 = l2norm_sqr(py - q.as_vec()) / hs2;
                      double d1 = l2norm_sqr(vy - f(q).as_vec()) / hr2;
                      double d = g(d0 + d1);
                      s1 += d * q.as_vec();
                      s2 += d * f(q).as_vec();
                      s3 += d;
                    }
                }
              site_t new_py = s1 / s3;
              value_t new_vy = s2 / s3;
              stop = (new_py == py);
              if (not stop)
                {
                  double d0 = l2norm_sqr(py - new_py);
                  double d1 = l2norm_sqr(vy - new_vy);
                  stop = (d0 + d1) < eps2;
                  //if (stop)
                  //  std::cout << "#NITER:" << i << p << new_py << std::endl;
                }
              py = new_py;
              vy = new_vy;
            }
          //std::cout << vy.as_vec() << std::endl;
          out(p) = (V)vy.as_vec();
        }
      return out;
    }


  } // end of namespace mln::contrib

} // end of namespace mln

#endif //!MLN_CONTRIB_MEANSHIFT_MEANSHIFT_HPP
