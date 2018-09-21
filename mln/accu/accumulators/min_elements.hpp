#ifndef MLN_ACCU_ACCUMULATORS_MIN_ELEMENTS_HPP
# define MLN_ACCU_ACCUMULATORS_MIN_ELEMENTS_HPP

/// \file
/// \brief Header file for the minimum accumulator

# include <mln/accu/accumulator.hpp>
# include <mln/accu/accumulator_base.hpp>
# include <mln/core/value/value_traits.hpp>
# include <utility>
# include <vector>

namespace mln
{

  namespace accu
  {
    namespace accumulators
    {

      /// \brief Minimal elements accumulator.
      /// \tparam T The argument.
      /// \tparam Compare The comparison function (non strict partial or total order).
      template <typename T, typename Compare = productorder_less_equal<T> >
      struct minimal_elements;
    }

    namespace features
    {
      /// \brief Minimal elements accumulator.
      /// \tparam T The argument type.
      /// \tparam Compare The comparison function (`void` results in productorder_less by default).
      template <typename Compare = void>
      struct minimal_elements;
    }

    namespace extractor
    {

      template <typename A>
      auto
      minimal_elements(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), features::minimal_elements<> ()) )
      {
        return extract(exact(acc), features::minimal_elements<> ());
      }

    }


    namespace features
    {
      template <typename Compare>
      struct minimal_elements : simple_feature< minimal_elements<Compare> >
      {
        minimal_elements(const Compare& cmp = Compare())
          : m_cmp(cmp)
        {
        }

        template <typename T>
        struct apply
        {
          typedef accumulators::minimal_elements<T, Compare> type;
        };

        template <typename T>
        accumulators::minimal_elements<T, Compare>
        make() const
        {
          return accumulators::minimal_elements<T, Compare> (m_cmp);
        }

      private: // dynamic parameter.
        Compare m_cmp;
      };

      namespace internal
      {
        template <typename T>
        using meta_minimal_elements = accumulators::minimal_elements<T>;
      }

      template <>
      struct minimal_elements<void>
        : simple_feature_facade< minimal_elements<void>, internal::meta_minimal_elements>
      {
      };
    }

    namespace accumulators
    {

      template <typename T, typename Compare>
      struct minimal_elements : accumulator_base< minimal_elements<T, Compare>, T, T, features::minimal_elements<> >
      {
        typedef T       argument_type;
        typedef std::vector<T> return_type;
        //typedef features::min<> feature;

        minimal_elements(const Compare& cmp = Compare())
          : m_cmp( cmp )
        {
        }

        void init()
        {
          m_val.clear();
        }

        void take(const T& v)
        {
          bool inserted = false;
          for (T& x: m_val)
            if (m_cmp(x, v))
              return;
            else if (m_cmp(v, x)) // new minimal element
              {
                x = v;
                inserted = true;
              }

          if (not inserted)
            m_val.push_back(v);
        }

        template <typename Other>
        void take(const Accumulator<Other>& other)
        {
          std::vector<T> vec = extractor::minimal_elements(other);
          for (T v: vec)
            this->take(v);
        }

        friend
        std::vector<T>
        extract(const minimal_elements& accu, features::minimal_elements<> )
        {
          return accu.m_val;
        }

      private:
        std::vector<T> m_val;
        Compare m_cmp;
      };

    }

  }

}

#endif // ! MLN_ACCU_ACCUMULATORS_MIN_ELEMENTS_HPP
