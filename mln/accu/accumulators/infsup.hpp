#ifndef MLN_ACCU_ACCUMULATORS_INFSUP_HPP
# define MLN_ACCU_ACCUMULATORS_INFSUP_HPP

# include <mln/accu/accumulator_base.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/ops.hpp>
# include <utility>

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {

      template <typename T, typename Compare = productorder_less<T> >
      struct infsup;

      template <typename T, typename Compare = productorder_less<T> >
      struct sup;

      template <typename T, typename Compare = productorder_less<T> >
      struct inf;
    }

    namespace features
    {
      template <typename Compare = void>
      struct inf;

      template <typename Compare = void>
      struct sup;
    }

    namespace extractor
    {

      template <typename A>
      auto
      inf(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::inf<>> ()));

      template <typename A>
      auto
      sup(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::sup<>> ()));

    }

    namespace features
    {

      template <typename Compare>
      struct inf : simple_feature< inf<Compare> >
      {
        inf(const Compare& cmp)
          : m_cmp(cmp)
        {
        }

        template <typename T>
        struct apply
        {
          typedef accumulators::inf<T, Compare> type;
        };

        template <typename T>
        accumulators::inf<T, Compare>
        make() const
        {
          return accumulators::inf<T, Compare>(m_cmp);
        }

      private:
        Compare m_cmp;
      };

      template <>
      struct inf<void> : simple_feature< inf<void> >
      {
        inf() = default;

        template <typename T>
        struct apply
        {
          typedef accumulators::inf<T> type;
        };

        template <typename T>
        accumulators::inf<T>
        make() const
        {
          return accumulators::inf<T>();
        }
      };

      template <typename Compare>
      struct sup : simple_feature< sup<Compare> >
      {
        sup(const Compare& cmp)
          : m_cmp(cmp)
        {
        }

        template <typename T>
        struct apply
        {
          typedef accumulators::sup<T, Compare> type;
        };

        template <typename T>
        accumulators::sup<T, Compare>
        make() const
        {
          return accumulators::sup<T, Compare>(m_cmp);
        }

      private:
        Compare m_cmp;
      };

      template <>
      struct sup<void> : simple_feature< sup<void> >
      {
        sup() = default;

        template <typename T>
        struct apply
        {
          typedef accumulators::sup<T> type;
        };

        template <typename T>
        accumulators::sup<T>
        make() const
        {
          return accumulators::sup<T>();
        }
      };

    }

    namespace extractor
    {

      template <typename A>
      auto
      inf(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::inf<>> ()) )
      {
        return extract(exact(acc), features::inf<> ());
      }

      template <typename A>
      auto
      sup(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::sup<>> ()) )
      {
        return extract(exact(acc), features::sup<> ());
      }

    }

    namespace accumulators
    {

      template <typename T, typename Compare>
      struct infsup : Accumulator< infsup<T, Compare> >
      {
        typedef T                               argument_type;
        typedef std::pair<T,T>			result_type;
        typedef boost::mpl::set<
          features::inf<>, features::sup<> >    provides;

        infsup(const Compare& cmp = Compare())
          : m_inf( value_traits<T, Compare>::sup() ),
            m_sup( value_traits<T, Compare>::inf() ),
            m_cmp( cmp )
        {
        }


        void init()
        {
          m_inf = value_traits<T, Compare>::sup();
          m_sup = value_traits<T, Compare>::inf();
        }

        void take(const T& v)
        {
          using mln::inf;
          using mln::sup;
          m_inf = inf(m_inf, v, m_cmp);
          m_sup = sup(m_sup, v, m_cmp);
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          using mln::inf;
          using mln::sup;
          T vinf = extractor::inf(other);
          T vsup = extractor::sup(other);
          m_inf = inf(m_inf, vinf, m_cmp);
          m_sup = sup(m_sup, vsup, m_cmp);
        }


        std::pair<T, T> to_result() const
        {
          return std::make_pair(m_inf, m_sup);
        }

        friend
        T extract(const infsup& accu, features::inf<> )
        {
          return accu.m_inf;
        }

        friend
        T extract(const infsup& accu, features::sup<> )
        {
          return accu.m_sup;
        }

      private:
        T	m_inf;
        T	m_sup;
        Compare m_cmp;
      };

      template <typename T, typename Compare>
      struct inf : accumulator_base< inf<T, Compare>, T, T, features::inf<> >
      {
        typedef T           argument_type;
        typedef T           result_type;

        inf(const Compare& cmp = Compare())
          : m_cmp( cmp ),
            m_inf( value_traits<T, Compare>::sup() )
        {
        }

        void init()
        {
          m_inf = value_traits<T, Compare>::sup();
        }

        void take(const T& v)
        {
          using mln::inf;
          m_inf = inf(m_inf, v, m_cmp);
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          using mln::inf;
          T vinf = extractor::inf(other);
          m_inf = inf(m_inf, vinf, m_cmp);
        }

        template <class dummy = bool>
        typename std::enable_if<std::is_same<T, bool>::value, dummy>::type
        stop() const
        {
          return m_inf == value_traits<bool, Compare>::inf();
        }


        friend
        T extract(const inf& accu, features::inf<> )
        {
          return accu.m_inf;
        }

      private:
        Compare m_cmp;
        T       m_inf;
      };

      template <typename T, typename Compare>
      struct sup : accumulator_base< sup<T, Compare>, T, T, features::sup<> >
      {
        typedef T           argument_type;
        typedef T			result_type;

        sup(const Compare& cmp = Compare())
          : m_cmp( cmp ),
            m_sup( value_traits<T, Compare>::inf() )
        {
        }

        void init()
        {
          m_sup = value_traits<T, Compare>::inf();
        }

        void take(const T& v)
        {
          using mln::sup;
          m_sup = sup(m_sup, v, m_cmp);
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          using mln::sup;
          T vsup = extractor::sup(other);
          m_sup = sup(m_sup, vsup, m_cmp);
        }

        template <class dummy = bool>
        typename std::enable_if<std::is_same<T, bool>::value, dummy>::type
        stop() const
        {
          return m_sup == value_traits<bool, Compare>::sup();
        }

        friend
        T extract(const sup& accu, features::sup<> )
        {
          return accu.m_sup;
        }

      private:
        Compare m_cmp;
        T       m_sup;
      };

    }

  }

}

#endif // ! MLN_ACCU_ACCUMULATORS_INFSUP_HPP
