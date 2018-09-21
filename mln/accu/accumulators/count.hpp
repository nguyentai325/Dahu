#ifndef MLN_ACCU_ACCUMULATOTS_COUNT_HPP
# define MLN_ACCU_ACCUMULATOTS_COUNT_HPP

/// \file
/// \brief Header file for counting accumulator.

# include <mln/accu/accumulator_base.hpp>
# include <mln/accu/dontcare.hpp>
# include <utility>

namespace mln
{

  namespace accu
  {
    namespace accumulators
    {
      /// \brief The counting accumulator.
      /// \tparam CountType is underlying the counter type (it must be incrementable,
      ///         and decrementable.
      template <typename CountType = std::size_t>
      struct count;
    }

    namespace features
    {
      template <typename CountType = std::size_t>
      struct count;
    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      count(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), std::declval<features::count<>> ()));

    }

    namespace features
    {
      template <typename CountType>
      struct count : simple_feature< count<CountType> >
      {
        template <typename T>
        struct apply
        {
          typedef accumulators::count<CountType> type;
        };

        template <typename T>
        accumulators::count<CountType>
        make() const
        {
          return accumulators::count<CountType> ();
        }
      };
    }

    namespace extractor
    {

      template <typename A>
      inline
      auto
      count(const Accumulator<A>& acc)
        -> decltype(extract(exact(acc), std::declval<features::count<>> ()))
      {
        return extract(exact(acc), features::count<> ());
      }

    }


    namespace accumulators
    {

      template <typename CountType>
      struct count : accumulator_base< count<CountType>, dontcare, CountType, features::count<> >
      {
        typedef dontcare                                argument_type;
        typedef CountType                               result_type;
        typedef boost::mpl::set< features::count<> >	provides;
        typedef std::true_type                          has_untake;

        count()
          : m_count (0)
        {
        }

        void init()
        {
          m_count = 0;
        }

        void take(const dontcare&)
        {
          ++m_count;
        }

        void untake(const dontcare&)
        {
          --m_count;
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          m_count += extractor::count(other);
        }

        friend
        result_type extract(const count& accu, features::count<> )
        {
          return accu.m_count;
        }

      private:
        CountType m_count;
      };

    }

  }

}

#endif // !MLN_ACCU_ACCUMULATOTS_COUNT_HPP
