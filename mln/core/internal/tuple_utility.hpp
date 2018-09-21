#ifndef TUPLE_UTILITY_HPP
# define TUPLE_UTILITY_HPP

# include <tuple>
# include <iostream>
# include <mln/core/internal/intseq.hpp>

namespace mln
{

  namespace internal
  {


    template < class F, class... TTypes>
    std::tuple< typename F::template apply<TTypes&>::type... >
    tuple_transform(std::tuple<TTypes...>& t, F f);

    template < class F, class... TTypes>
    std::tuple< typename F::template apply<const TTypes&>::type... >
    tuple_transform(const std::tuple<TTypes...>& t, F f);

    template < class F, class... TTypes>
    std::tuple<TTypes...>&
    tuple_for_each(std::tuple<TTypes...>& t, F f);

    template < class F, class... TTypes>
    const std::tuple<TTypes...>&
    tuple_for_each(const std::tuple<TTypes...>& t, F f);

    template <class Tuple, class F>
    struct tuple_meta_transform;


    /************************/
    /*** Implementation   ***/
    /************************/

    namespace impl
    {

      // template <class THead, class... TTypes>
      // using head_t = THead;

      // template <class THead, class... TTail>
      // using tail_t = std::tuple<TTail...>;

      // template <int N, class F, class... TTypes>
      // struct tuple_transform_helper;

      // template <int N, class F, class Head, class... TTypes>
      // struct tuple_transform_helper<N, F, Head, TTypes...> : tuple_transform_helper<N-1, F, TTypes...>
      // {
      // };

      // template <class F, class Head, class... TTypes>
      // struct tuple_transform_helper<0, F, Head, TTypes...>
      // {
      //   typedef std::tuple< typename F::template apply<Head&>::type,
      //   		    typename F::template apply<TTypes&>::type... > type;
      // };

      // template <class F>
      // struct tuple_transform_helper<0, F>
      // {
      //   typedef std::tuple<> type;
      // };

      // template <int N, class F, class... TTypes>
      // struct const_tuple_transform_helper;

      // template <int N, class F, class Head, class... TTypes>
      // struct const_tuple_transform_helper<N, F, Head, TTypes...> :
      //   const_tuple_transform_helper<N-1, F, TTypes...>
      // {
      // };

      // template <class F, class Head, class... TTypes>
      // struct const_tuple_transform_helper<0, F, Head, TTypes...>
      // {
      //   typedef std::tuple< typename F::template apply<const Head&>::type,
      //   		    typename F::template apply<const TTypes&>::type... > type;
      // };

      // template <class F>
      // struct const_tuple_transform_helper<0, F>
      // {
      //   typedef std::tuple<> type;
      // };


      // Transform: terminal case
      // template <int N, class F, class... TTypes>
      // std::tuple<>
      // tuple_transform(const std::tuple<TTypes...>&, F,
      //   	      typename std::enable_if<N == sizeof...(TTypes)>::type* = 0)
      // {
      //   return std::tuple<>();
      // }


      template <class F, class...>
      struct tuple_transform_helper;

      template <class F, class...>
      struct const_tuple_transform_helper;


      template <class F, int... I, class... TTypes>
      struct tuple_transform_helper<F, intseq<I...>, TTypes...>
      {
        typedef std::tuple<
          typename std::result_of<F(TTypes&)>::type...
          > result_type;

        result_type
        operator() (std::tuple<TTypes...>& t, F f)
        {
          return result_type(f(std::get<I>(t))...);
        }
      };


      template <class F, int... I, class... TTypes>
      struct const_tuple_transform_helper<F, intseq<I...>, TTypes...>
      {
        typedef std::tuple<
          typename std::result_of<F(const TTypes&)>::type...
          > result_type;

        result_type
        operator() (const std::tuple<TTypes...>& t, F f)
        {
          return result_type(f(std::get<I>(t))...);
        }
      };




      // // Transform: rec case (non-const)
      // template <int N, class F, class... TTypes>
      // typename tuple_transform_helper<N, F, TTypes...>::type
      // tuple_transform(std::tuple<TTypes...>& t, F f,
      //   	      typename std::enable_if<N != sizeof...(TTypes)>::type* = 0)
      // {
      //   typedef typename std::tuple_element<N, std::tuple<TTypes...> >::type _head;
      //   typedef typename F::template apply<_head&>::type THead;

      //   std::tuple<THead> head(f(std::get<N>(t)));
      //   auto tail = tuple_transform<N+1>(t, f);
      //   auto x = std::tuple_cat(head, tail);
      //   return x;
      // }

      // // Transform: rec case (const)
      // template <int N, class F, class... TTypes>
      // typename const_tuple_transform_helper<N, F, TTypes...>::type
      // tuple_transform(const std::tuple<TTypes...>& t, F f,
      //   	      typename std::enable_if<N != sizeof...(TTypes)>::type* = 0)
      // {
      //   typedef typename std::tuple_element<N, std::tuple<TTypes...> >::type _head;
      //   typedef typename F::template apply<const _head&>::type THead;

      //   std::tuple<THead> head(f(std::get<N>(t)));
      //   auto tail = tuple_transform<N+1>(t, f);
      //   auto x = std::tuple_cat(head, tail);
      //   return x;
      // }

      // For each: terminal case
      template <int N, class F, class... TTypes>
      void
      tuple_for_each(const std::tuple<TTypes...>&, F,
		     typename std::enable_if<N == sizeof...(TTypes)>::type* = 0)
      {
      }

      // For each: rec case (non-const)
      template <int N, class F, class... TTypes>
      void
      tuple_for_each(std::tuple<TTypes...>& t, F f,
		     typename std::enable_if<N != sizeof...(TTypes)>::type* = 0)
      {
	f(std::get<N>(t));
	tuple_for_each<N+1>(t, f);
      }

      // For each: rec case (const)
      template <int N, class F, class... TTypes>
      void
      tuple_for_each(const std::tuple<TTypes...>& t, F f,
		     typename std::enable_if<N != sizeof...(TTypes)>::type* = 0)
      {
	f(std::get<N>(t));
	tuple_for_each<N+1>(t, f);
      }

    }

    template <class F, class... TTypes>
    std::tuple< typename F::template apply<TTypes&>::type... >
    tuple_transform(std::tuple<TTypes...>& t, F f)
    {
      typedef typename int_list_seq<sizeof...(TTypes)>::type S;
      return impl::tuple_transform_helper<F, S, TTypes...> () (t, f);
    }

    template <class F, class... TTypes>
    std::tuple< typename F::template apply<const TTypes&>::type... >
    tuple_transform(const std::tuple<TTypes...>& t, F f)
    {
      typedef typename int_list_seq<sizeof...(TTypes)>::type S;
      return impl::const_tuple_transform_helper<F, S, TTypes...> () (t, f);
    }



    template <class F, class... TTypes>
    std::tuple<TTypes...>&
    tuple_for_each(std::tuple<TTypes...>& t, F f)
    {
      impl::tuple_for_each<0>(t, f);
      return t;
    }

    template <class F, class... TTypes>
    const std::tuple<TTypes...>&
    tuple_for_each(const std::tuple<TTypes...>& t, F f)
    {
      impl::tuple_for_each<0>(t, f);
      return t;
    }


    template <class F, class... TTypes>
    struct tuple_meta_transform< std::tuple<TTypes...>, F>
    {
      typedef std::tuple< typename F::template apply<TTypes&>::type... > type;
    };

    template <class F, class... TTypes>
    struct tuple_meta_transform< const std::tuple<TTypes...>, F>
    {
      typedef std::tuple< typename F::template apply<const TTypes&>::type... > type;
    };

  }

}

#endif // ! TUPLE_UTILITY_HPP
