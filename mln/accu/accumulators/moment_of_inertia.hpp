#ifndef MLN_ACCU_ACCUMULATORS_MOMENT_OF_INERTIA_HPP
# define MLN_ACCU_ACCUMULATORS_MOMENT_OF_INERTIA_HPP

# include <mln/core/math_ops.hpp>
# include <mln/accu/accumulators/sum.hpp>
# include <mln/accu/accumulators/count.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/accu/composite_accumulator.hpp>
# include <boost/type_traits/promote.hpp>

/// \file
/// \brief Moment of inertia accumulator.

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {

      /// \brief Accumulator for the scale invariant moment of inertia
      ///
      /// The moment of inertia is the 1st hu moment, i.e., in the 2D
      /// case:
      ///
      /// \f$ \frac{1}{n^(1+2/d)} \mu_{(2,0)} + \mu_{(0,2)} \f$
      /// with:
      /// \f[
      ///  \mu_{pq} = \sum_{x} \sum_{y} (x - \bar{x})^p(y -
      /// \bar{y})^q
      /// \f]
      template <typename T,
                typename SumType = typename boost::promote<T>::type,
                typename SumSqrType = SumType>
      struct moment_of_inertia;
    }

    namespace features
    {
      template <typename SumType = void,
                typename SumSqrType = SumType>
      struct moment_of_inertia;
    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      moment_of_inertia(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), features::moment_of_inertia<> ()))
      {
        return extract(exact(acc), features::moment_of_inertia<> ());
      }

    }


    namespace features
    {
      namespace internal
      {
        template <typename SumType, typename SumSqrType>
        struct meta_moment_of_inertia
        {
          template <typename T>
          using type = accu::accumulators::moment_of_inertia<T, SumType, SumSqrType>;
        };

        template <>
        struct meta_moment_of_inertia<void, void>
        {
          template <typename T>
          using type = accu::accumulators::moment_of_inertia<T>;
        };
      }

      template <typename SumType,
                typename SumSqrType>
      struct moment_of_inertia : simple_feature_facade<
        moment_of_inertia<SumType, SumSqrType>,
        internal::meta_moment_of_inertia<SumType, SumSqrType>::template type
        >
      {
      };
    }

    namespace accumulators
    {

      template <typename T, typename SumType, typename SumSqrType>
      struct moment_of_inertia : accumulator_base< moment_of_inertia<T, SumType, SumSqrType>,
                                          T,
                                          double,
                                          features::moment_of_inertia<> >
      {
        typedef T       argument_type;
        typedef double  result_type;


        moment_of_inertia()
          : m_count {0},
            m_sum {},
            m_sum_sqr {}
        {
        }

        void init()
        {
          m_count = 0;
          m_sum = SumType ();
          m_sum_sqr = SumSqrType ();
        }

        void take(const argument_type& arg)
        {
          SumSqrType x(arg);
          m_count += 1;
          m_sum += SumType(arg);
          m_sum_sqr += x * x;
        }

        void untake(const argument_type& arg)
        {
          SumSqrType x(arg);
          m_count -= 1;
          m_sum -= SumType(arg);
          m_sum_sqr -= x * x;
        }

        void take(const moment_of_inertia& other)
        {
          m_count += other.m_count;
          m_sum += other.m_sum;
          m_sum_sqr += other.m_sum_sqr;
        }

        void untake(const moment_of_inertia& other)
        {
          m_count -= other.m_count;
          m_sum -= other.m_sum;
          m_sum_sqr -= other.m_sum_sqr;
        }


        template <class ST1, class ST2>
        friend
        result_type extract(const moment_of_inertia& accu, features::moment_of_inertia<ST1, ST2> )
        {
          using mln::sum;

          if (accu.m_count == 0)
            return 0;

          result_type n = accu.m_count;
          result_type inertia = sum(accu.m_sum_sqr - sqr(accu.m_sum) / n);

          constexpr unsigned dim = value_traits<typename moment_of_inertia::argument_type>::ndim;
          if (dim == 2)
            return inertia / sqr(n);
          else
            return inertia / std::pow(n, 1.0 + 2.0 / value_traits<typename moment_of_inertia::argument_type>::ndim);
        }

        friend
        unsigned extract(const moment_of_inertia& accu, features::count<> )
        {
          return accu.m_count;
        }

        friend
        SumType
        extract(const moment_of_inertia& accu, features::mean<> )
        {
          return accu.m_sum / accu.m_count;
        }

      private:
        unsigned        m_count;
        SumType         m_sum;
        SumSqrType      m_sum_sqr;
      };

    }

  }

}

#endif // ! MLN_ACCU_ACCUMULATORS_MOMENT_OF_INERTIA_HPP
