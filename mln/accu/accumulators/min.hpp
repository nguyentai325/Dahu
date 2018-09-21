#ifndef MLN_ACCU_ACCUMULATORS_MIN_HPP
# define MLN_ACCU_ACCUMULATORS_MIN_HPP

/// \file
/// FIXME: use literal::zero instead of default initialization

# include <mln/accu/accumulator_base.hpp>
# include <mln/core/value/value_traits.hpp>
# include <utility>

namespace mln
{

  namespace accu
  {
    namespace accumulators
    {

      template <typename T, typename Compare = std::less<T> >
      struct min;
    }

    namespace features
    {
      template <typename Compare = void>
      struct min;
    }

    namespace extractor
    {

      template <typename A>
      auto
      min(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::min<>> ()) );

    }


    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace features
    {

      template <typename Compare>
      struct min : simple_feature< min<Compare> >
      {
        min(const Compare& cmp)
          : m_cmp(cmp)
        {
        }

        template <typename T>
        struct apply
        {
          typedef accumulators::min<T, Compare> type;
        };

        template <typename T>
        accumulators::min<T, Compare>
        make() const
        {
          return accumulators::min<T, Compare>(m_cmp);
        }

      private:
        Compare m_cmp;
      };

      namespace internal
      {
        template <typename T>
        using meta_min = accumulators::min<T>;
      }

      template <>
      struct min<void>
        : simple_feature_facade< min<void>, internal::meta_min>
      {
      };
    }

    namespace extractor
    {

      template <typename A>
      auto
      min(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), std::declval<features::min<>> ()) )
      {
        return extract(exact(acc), features::min<> ());
      }

    }

    namespace accumulators
    {

      template <typename T, typename Compare>
      struct min : accumulator_base< min<T, Compare>, T, T, features::min<> >
      {
        typedef T       argument_type;
        typedef T       return_type;
        //typedef features::min<> feature;

        min(const Compare& cmp = Compare())
          : m_val( value_traits<T, Compare>::max() ),
            m_cmp( cmp )
        {
        }

        void init()
        {
          m_val = value_traits<T, Compare>::max();
        }

        void take(const T& v)
        {
          if (m_cmp(v, m_val))
            m_val = v;
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          T v = extractor::min(other);
          if (m_cmp(v, m_val))
            m_val = v;
        }

        friend
        T extract(const min& accu, features::min<> )
        {
          return accu.m_val;
        }

      private:
        T   m_val;
        Compare m_cmp;
      };

    }

  }

}

#endif // !MLN_ACCU_ACCUMULATORS_MIN_HPP
