#ifndef MLN_ACCU_ACCUMULATORS_ACCU_AS_IT_HPP
# define MLN_ACCU_ACCUMULATORS_ACCU_AS_IT_HPP

# include <mln/accu/accumulator_base.hpp>

namespace mln
{

  namespace accu
  {

    namespace accumulators
    {

      /// \brief Wrapper around an accumlator that forces the
      /// result type to be the accu itself instead of a result_type
      template <class Accu>
      struct accu_as_it
        : Accumulator< accu_as_it<Accu> >
      {
        typedef typename Accu::argument_type argument_type;
        typedef Accu                         result_type;

        accu_as_it() = default;

        accu_as_it(const Accu& x)
          : m_accu(x)
        {
        }

        accu_as_it(Accu&& x)
          : m_accu(std::move(x))
        {
        }

        void init()
        {
          m_accu.init();
        }

        void take(const argument_type& x)
        {
          m_accu.take(x);
        }

        void take(const Accu& x)
        {
          m_accu.take(x);
        }

        void take(const accu_as_it& x)
        {
          m_accu.take(x.m_accu);
        }

        void untake(const argument_type& x)
        {
          m_accu.untake(x);
        }


        void untake(const Accu& x)
        {
          m_accu.untake(x.m_accu);
        }

        const Accu& to_result() const
        {
          return m_accu;
        }

      private:
        Accu m_accu;
      };
    }
  }
}

#endif // ! MLN_ACCU_ACCUMULATORS_ACCU_AS_IT_HPP
