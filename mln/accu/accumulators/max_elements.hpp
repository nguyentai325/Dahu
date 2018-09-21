#ifndef MLN_ACCU_ACCUMULATORS_MAX_ELEMENTS_HPP
# define MLN_ACCU_ACCUMULATORS_MAX_ELEMENTS_HPP

/// \file
/// \brief Header file for the maximum accumulator

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

      /// \brief Maximal elements accumulator.
      ///
      /// Given a set (V,≺) where ≺ denotes a partial order.
      /// x is a maximal element of a subset S ⊆ V if there is no
      /// other y ∈ S that is greater than x. More formally:
      /// ∀ y ∈ S, ¬(x ≺ y).
      ///
      /// Note that the returned set may contain doublons.
      /// \tparam T The type of elements to accumulate
      /// \tparam Compare The comparison function ≼ (non-strict).
      template <typename T, typename Compare = productorder_less_equal<T> >
      struct maximal_elements;
    }

    namespace features
    {
      /// \brief Maximal elements accumulator.
      /// \tparam T The argument type.
      /// \tparam Compare The comparison function (`void` results in productorder_less by default).
      template <typename Compare = void>
      struct maximal_elements;
    }

    namespace extractor
    {

      template <typename A>
      auto
      maximal_elements(const Accumulator<A>& acc)
        -> decltype( extract(exact(acc), features::maximal_elements<> ()) )
      {
        return extract(exact(acc), features::maximal_elements<> ());
      }

    }


    namespace features
    {
      template <typename Compare>
      struct maximal_elements : simple_feature< maximal_elements<Compare> >
      {
        maximal_elements(const Compare& cmp = Compare())
          : m_cmp(cmp)
        {
        }

        template <typename T>
        struct apply
        {
          typedef accumulators::maximal_elements<T, Compare> type;
        };

        template <typename T>
        accumulators::maximal_elements<T, Compare>
        make() const
        {
          return accumulators::maximal_elements<T, Compare> (m_cmp);
        }

      private: // dynamic parameter.
        Compare m_cmp;
      };

      namespace internal
      {
        template <typename T>
        using meta_maximal_elements = accumulators::maximal_elements<T>;
      }

      template <>
      struct maximal_elements<void>
        : simple_feature_facade< maximal_elements<void>, internal::meta_maximal_elements>
      {
      };
    }

    namespace accumulators
    {

      template <typename T, typename Compare>
      struct maximal_elements : accumulator_base< maximal_elements<T, Compare>, T, T,
                                                  features::maximal_elements<> >
      {
        typedef T       argument_type;
        typedef std::vector<T> return_type;
        //typedef features::max<> feature;

        maximal_elements(const Compare& cmp = Compare())
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
            if (m_cmp(v, x))
              return;
            else if (m_cmp(x, v)) // new maximal element
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
          std::vector<T> vec = extractor::maximal_elements(other);
          for (T v: vec)
            this->take(v);
        }

        friend
        std::vector<T>
        extract(const maximal_elements& accu, features::maximal_elements<> )
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

#endif // ! MLN_ACCU_ACCUMULATORS_MAX_ELEMENTS_HPP
