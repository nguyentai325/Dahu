#ifndef COMPOSITE_ACCUMULATOR_HPP
# define COMPOSITE_ACCUMULATOR_HPP

# include <mln/accu/concept/accumulator.hpp>
# include <mln/accu/feature.hpp>
# include <boost/mpl/assert.hpp>
# include <boost/mpl/set.hpp>
//# include <boost/mpl/set10.hpp>
# include <boost/mpl/transform.hpp>
# include <boost/mpl/fold.hpp>
# include <boost/mpl/copy.hpp>
# include <boost/mpl/inserter.hpp>
# include <boost/mpl/remove_if.hpp>
# include <boost/mpl/placeholders.hpp>
# include <boost/mpl/quote.hpp>
# include <boost/mpl/insert.hpp>
# include <boost/mpl/has_key.hpp>
# include <boost/mpl/empty.hpp>
# include <boost/mpl/protect.hpp>
# include <boost/fusion/mpl.hpp>

namespace mln
{

  namespace accu
  {

    namespace internal
    {

      namespace mpl = boost::mpl;
      namespace ph = boost::mpl::placeholders;

      template <typename Set>
      using mpl_set_inserter = mpl::inserter< Set, mpl::insert<ph::_1, ph::_2> >;

      template <typename Set>
      using mpl_set_inserter_p = mpl::inserter< Set, mpl::lambda< mpl::insert<ph::_, ph::_> > >;

      template <typename Set, typename Entry>
      struct whatdepends_helper_copy
        : mpl::copy<Entry, mpl_set_inserter<Set> >
      {
      };

      template <typename FeatureSet>
      struct whatdepends
      {

        typedef mpl::quote1<features::depends>     meta_depends;
        typedef typename mpl::transform<FeatureSet, meta_depends, mpl_set_inserter< mpl::set<> > >::type dependances;
        typedef typename mpl::fold< dependances, mpl::set<>,
                                    whatdepends_helper_copy<ph::_1, ph::_2>
                                    >::type type;
      };

      template <typename FeatureSet, typename Depends>
      struct get_unresolved
      {
        typedef typename mpl::remove_if< Depends, mpl::has_key<FeatureSet, ph::_1>, mpl_set_inserter< mpl::set<> > >::type type;
      };


      template <typename FeatureSet, typename Depends = typename whatdepends<FeatureSet>::type,
                bool empty = mpl::empty<Depends>::value >
      struct resolve_dependances
      {
        typedef typename get_unresolved<FeatureSet, Depends>::type unresolved;
        typedef typename mpl::copy<unresolved, mpl_set_inserter<FeatureSet> >::type new_feature_set;
        typedef typename resolve_dependances<new_feature_set, typename whatdepends<unresolved>::type>::type type;
      };

      template <typename FeatureSet, typename Depends>
      struct resolve_dependances<FeatureSet, Depends, true>
      {
        typedef FeatureSet type;
      };


      template <typename T>
      struct feature_to_accu_helper
      {
        template <typename F>
        struct apply
        {
          typedef typename F::template apply<T>::type type;
        };
      };

      struct accu_init
      {
        template <typename Accu>
        void operator () (Accu& acc) const { acc.init(); }
      };

      template <typename T>
      struct accu_take
      {
        accu_take(const T& x) : m_x(x) {}

        template <typename Accu>
        void operator () (Accu& acc) const { acc.take(m_x); }

      private:
        const T& m_x;
      };

      template <typename T>
      struct accu_untake
      {
        accu_untake(const T& x) : m_x(x) {}

        template <typename Accu>
        void operator () (Accu& acc) const { acc.untake(m_x); }

      private:
        const T& m_x;
      };

      template <typename Feature>
      struct accu_has_feature
      {
        template <typename Accu>
        struct apply
        {
          typedef typename boost::mpl::has_key<typename Accu::provides, Feature>::type type;
        };
      };


      template <typename AccuList, typename Feature>
      struct acculist_has_feature
      {
        typedef typename boost::fusion::result_of::find_if<AccuList, accu_has_feature<Feature> >::type iterator;
        typedef typename boost::fusion::result_of::end<AccuList>::type end;
        static constexpr bool value = not boost::fusion::result_of::equal_to<iterator, end>::value;


        //static constexpr iterator x = end ();
      };

      template <typename AccuList, typename Feature>
      struct acculist_get_feature
      {
        typedef typename boost::fusion::result_of::find_if<AccuList, accu_has_feature<Feature> >::type iterator;
        typedef typename boost::fusion::result_of::deref<iterator>::type accu;

	
        typedef decltype(extract (std::declval<accu>(), Feature ()) ) type;
      };
    }


    template <typename E, typename T, typename FeatureSet>
    struct composite_accumulator_base : Accumulator<E>
    {
      typedef T					 argument_type;
      typedef typename internal::resolve_dependances<typename FeatureSet::features>::type	        fset;
      typedef typename boost::mpl::transform<fset, internal::feature_to_accu_helper<T>,
                                             boost::mpl::back_inserter< boost::fusion::list<> > >::type acculist;

      typedef E result_type;

      void init()
      {
        boost::fusion::for_each(m_accus, internal::accu_init ());
      }


      void take(const argument_type& x)
      {
        boost::fusion::for_each(m_accus, internal::accu_take<argument_type>(x));
      }

      template <typename Other>
      void take(const Accumulator<Other>& other)
      {
        boost::fusion::for_each(m_accus, internal::accu_take<Other>(exact(other)));
      }


      void untake(const argument_type& x)
      {
        boost::fusion::for_each(m_accus, internal::accu_untake<argument_type>(x));
      }

      template <typename Feature>
      friend
      typename boost::lazy_enable_if_c< internal::acculist_has_feature<acculist, Feature>::value,
                                        internal::acculist_get_feature<acculist, Feature> >::type
      extract(const composite_accumulator_base& accu, Feature feat)
      {
        auto res = boost::fusion::find_if< internal::accu_has_feature<Feature> >(accu.m_accus);

        return extract(*res, feat);
      }

      const E& to_result() const
      {
        return exact(*this);
      }

    private:
      acculist m_accus;
    };

    template <typename T, typename FeatureSet>
    struct composite_accumulator :
      composite_accumulator_base< composite_accumulator<T, FeatureSet>, T, FeatureSet>
    {
    };

    template <typename E, typename ArgumentType, typename ResultType, typename Feature>
    struct composite_accumulator_facade :
      composite_accumulator_base<E, ArgumentType, features::composite_feature<typename features::depends<Feature>::type> >
    {
      typedef ArgumentType	   argument_type;
      typedef ResultType	   result_type;
      typedef Feature		   feature;

      ResultType to_result() const
      {
        return extract(exact(*this), Feature() );
      }
    };

  }

}

#endif // ! COMPOSITE_ACCUMULATOR_HPP
