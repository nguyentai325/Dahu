#ifndef MLN_ACCU_ACCUMULATORS_H_INFSUP_HPP
# define MLN_ACCU_ACCUMULATORS_H_INFSUP_HPP

# include <mln/accu/accumulator_base.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/value/indexer.hpp>
# include <mln/accu/accumulators/infsup.hpp>

/// FIXME: use indexer / rename has histogram

namespace mln
{

  namespace accu
  {


    namespace accumulators
    {
      /// \brief An accumulator that computes the infimum of values
      /// based on the tracking of values. It provides 
      template <class T>
      struct h_inf;

      template <class T>
      struct h_sup;
    }

    namespace features
    {
      struct h_inf;
      struct h_sup;
    }

    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace features
    {

      struct h_inf : simple_feature_facade<h_inf, accumulators::h_inf>
      {
      };

      struct h_sup : simple_feature_facade<h_sup, accumulators::h_sup>
      {
      };

    }


    namespace accumulators
    {

      template <typename E, typename T, typename F, typename Enable = void>
      struct h_infsup_base;


      template <typename E, typename T, typename F>
      struct h_infsup_base<E, T, F, typename std::enable_if<
                                      std::is_integral<T>::value and
                                      value_traits<T>::quant <= 16>::type>
      : accumulator_base<E, T, T, F>
      {
        typedef T argument_type;
        typedef T result_type;
        typedef boost::mpl::set<
          features::h_sup,
          features::h_inf,
          features::inf<>,
          features::sup<> > provides;

        h_infsup_base()
          : m_inf(value_traits<T>::sup()),
            m_sup(value_traits<T>::inf()),
            m_count (0),
            m_hist {{0,}}
        {
        }

        void init()
        {
          m_inf = value_traits<T>::sup();
          m_sup = value_traits<T>::inf();
          m_count = 0;
          m_hist.fill(0);
        }

        void take(const T& x)
        {
          ++m_hist[x];
          ++m_count;
          if (x < m_inf) m_inf = x;
          if (x > m_sup) m_sup = x;
        }

        void untake(const T& x)
        {
          mln_precondition(m_hist[x] > 0);
          --m_hist[x];
          --m_count;

          if (m_hist[x] == 0)
            {
              if (m_count == 0) {
                m_inf = value_traits<T>::sup();
                m_sup = value_traits<T>::inf();
              } else if (x == m_inf) {
                int i = m_inf;
                while (not m_hist[i])
                  ++i;
                m_inf = i;
              } else if (x == m_sup) {
                int i = m_sup;
                while (not m_hist[i])
                  --i;
                m_sup = i;
              }
            }
        }

        friend
        T extract(const h_infsup_base& accu, features::inf<>)
        {
          return accu.m_inf;
        }

        friend
        T extract(const h_infsup_base& accu, features::sup<>)
        {
          return accu.m_sup;
        }

      protected:
        T m_inf;
        T m_sup;
        unsigned m_count;
        std::array<unsigned, indexer<T, std::less<T>>::nvalues> m_hist;
      };

      template <typename T>
      struct h_inf : h_infsup_base<h_inf<T>, T, features::inf<> >
      {
      };

      template <typename T>
      struct h_sup : h_infsup_base<h_sup<T>, T, features::sup<> >
      {
      };

    }
  }
}


#endif //!MLN_ACCU_ACCUMULATORS_H_INFSUP_HPP
