#ifndef MLN_ACCU_ACCUMULATORS_VARIANCE_HPP
# define MLN_ACCU_ACCUMULATORS_VARIANCE_HPP

# include <mln/accu/accumulators/sum.hpp>
# include <mln/accu/accumulators/count.hpp>
# include <mln/accu/accumulators/mean.hpp>
# include <mln/accu/composite_accumulator.hpp>
# include <boost/type_traits/promote.hpp>

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {

      /// \brief Accumulator for the variance
      /// \f[
      /// V(X) = E(\norm X - E(X) \norm_2^2)
      /// \f]
      template <typename T,
                typename SumType = typename boost::promote<T>::type,
                typename SumSqrType = SumType>
      struct variance;
    }

    namespace features
    {
      template <typename SumType = void,
                typename SumSqrType = SumType>
      struct variance;
    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      variance(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), features::variance<> ()))
      {
        return extract(exact(acc), features::variance<> ());
      }

    }


    namespace features
    {
      namespace internal
      {
        template <typename SumType, typename SumSqrType>
        struct meta_variance
        {
          template <typename T>
          using type = accu::accumulators::variance<T, SumType, SumSqrType>;
        };

        template <>
        struct meta_variance<void, void>
        {
          template <typename T>
          using type = accu::accumulators::variance<T>;
        };
      }

      template <typename SumType,
                typename SumSqrType>
      struct variance : simple_feature_facade< variance<SumType, SumSqrType>,
                                               internal::meta_variance<SumType, SumSqrType>::template type>
      {
      };
    }

    namespace accumulators
    {

      template <typename T, typename SumType, typename SumSqrType>
      struct variance : accumulator_base< variance<T, SumType, SumSqrType>,
                                          T,
                                          double,
                                          features::variance<> >
      {
        typedef T       argument_type;
        //typedef typename std::common_type<SumType, SumSqrType>::type result_type;
        typedef double  result_type;

        //typedef boost::mpl::set< features::variance<>, features::variance<SumType> > provides;

        variance()
          : m_count {0},
          m_sum {},
          m_sum_sqr {}
        {
        }

        void init()
        {
          m_count = 0;
          m_sum = SumType ();
          m_sum_sqr = result_type ();
        }

        void take(const argument_type& arg)
        {
          m_count += 1;
          m_sum += SumType(arg);
          m_sum_sqr += l2norm_sqr(SumSqrType(arg));
        }

        void untake(const argument_type& arg)
        {
          m_count -= 1;
          m_sum -= SumType(arg);
          m_sum_sqr -= l2norm_sqr(SumSqrType(arg));
        }

        void take(const variance& other)
        {
          m_count += other.m_count;
          m_sum += other.m_sum;
          m_sum_sqr += other.m_sum_sqr;
        }

        void untake(const variance& other)
        {
          m_count -= other.m_count;
          m_sum -= other.m_sum;
          m_sum_sqr -= other.m_sum_sqr;
        }


        template <class ST1, class ST2>
        friend
        result_type extract(const variance& accu, features::variance<ST1, ST2> )
        {
          if (accu.m_count == 0)
            return 0;

          result_type n = accu.m_count;
          return (accu.m_sum_sqr / n) - l2norm_sqr(accu.m_sum / n);
        }

        friend
        unsigned extract(const variance& accu, features::count<> )
        {
          return accu.m_count;
        }

        friend
        SumType
        extract(const variance& accu, features::mean<> )
        {
          return accu.m_sum / accu.m_count;
        }

      private:
        unsigned        m_count;
        SumType         m_sum;
        result_type     m_sum_sqr;
      };

    }

  }

}

#endif // ! MLN_ACCU_ACCUMULATORS_VARIANCE_HPP_HPP
