#ifndef MLN_ACCU_ACCUMULATOTS_MEAN_HPP
# define MLN_ACCU_ACCUMULATOTS_MEAN_HPP

# include <mln/accu/accumulators/sum.hpp>
# include <mln/accu/accumulators/count.hpp>
# include <mln/accu/composite_accumulator.hpp>

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {
      template <typename T, typename SumType = decltype( std::declval<T>() + std::declval<T>() )  >
      struct mean;
    }

    namespace features
    {
      template <typename SumType = void>
      struct mean;
    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      mean(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), std::declval<features::mean<>> ()));

    }


    namespace features
    {
      template <typename SumType>
      struct mean : simple_feature< mean<SumType> >
      {
        template <typename T>
        struct apply
        {
          typedef accumulators::mean<T, SumType> type;
        };

        template <typename T>
        accumulators::mean<T, SumType>
        make() const
        {
          return accumulators::mean<T, SumType>();
        }
      };

      template <>
      struct mean<void> : simple_feature< mean<void> >
      {
        template <typename T>
        struct apply
        {
          typedef accumulators::mean<T> type;
        };

        template <typename T>
        accumulators::mean<T>
        make() const
        {
          return accumulators::mean<T>();
        }
      };

      template <typename SumType>
      struct depends< mean<SumType> >
      {
        typedef boost::mpl::set< sum<SumType>, count<> > type;
      };

    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      mean(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), std::declval<features::mean<>> ()))
      {
        return extract(exact(acc), features::mean<> ());
      }

    }


    namespace accumulators
    {

      template <typename T, typename SumType>
      struct mean : composite_accumulator_facade< mean<T, SumType>, T,  SumType, features::mean<SumType> >
      {
        typedef T       argument_type;
        typedef SumType result_type;
        typedef boost::mpl::set< features::mean<>, features::mean<SumType> > provides;

        friend
        SumType extract(const mean& accu, features::mean<SumType> )
        {
          return extract(accu, features::mean<>() );
        }

        friend
        SumType extract(const mean& accu, features::mean<> )
        {
          auto v = extractor::count(accu);
          if (v == 0)
            return extractor::sum(accu);
          else
            return extractor::sum(accu) / extractor::count(accu);
        }
      };

    }

  }

}

#endif // ! MLN_ACCU_ACCUMULATOTS_MEAN_HPP
