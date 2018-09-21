#ifndef MLN_KERNELV2_DETAILS_EXPRESSIONS_TRAITS_HPP
# define MLN_KERNELV2_DETAILS_EXPRESSIONS_TRAITS_HPP

# include <boost/proto/proto.hpp>
# include <boost/proto/fusion.hpp>
# include <boost/mpl/list.hpp>
# include <boost/mpl/max_element.hpp>
# include <boost/mpl/fold.hpp>
# include <boost/mpl/int.hpp>
# include <boost/mpl/bool.hpp>
# include <boost/mpl/or.hpp>
# include <boost/mpl/if.hpp>
# include <boost/mpl/greater.hpp>
# include <boost/mpl/void.hpp>
# include <boost/type_traits.hpp>
# include <boost/fusion/tuple.hpp>
# include <mln/core/internal/intseq.hpp>
# include <mln/core/dontcare.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/image/constant_image.hpp>

namespace mln
{
  namespace kernel
  {
    namespace details
    {
      namespace proto = boost::proto;
      namespace mpl = boost::mpl;

      /// \brief Get expression traits
      /// It defines:
      /// * number_of_image (mpl::int)
      template <class Expr,
                class Tag = typename proto::tag_of<Expr>::type>
      struct expression_traits;

      template <class Expr>
      using number_of_images = typename expression_traits<Expr>::number_of_images;

      /// \brief Traits to get image usage
      /// It defines:
      /// * type The type of the image
      /// * point (mpl::bool_): true if f(p) exists
      /// * neighbor (mpl::bool_): true if f(n) exists
      template <class Expr, int i,
                class Tag = typename proto::tag_of<Expr>::type>
      struct image_usage_traits;


      /// \brief Get the image tuple from an expression
      template <class Expr>
      struct image_list_traits;

      template <class Expr>
      struct image_used_by_neighbor_list_traits;

      template<class Expr>
      typename image_list_traits<Expr>::type
      get_image_list(Expr&& x);

      template<class Expr>
      typename image_list_traits<Expr>::zip_image_type
      get_image_list2(Expr&& x);

      /// \brief Return a tuple of images from the expression
      /// where image that were not used as F(n) are replaced
      /// by a Null Image.
      /// A null image is of type constant_image<Domain, dontcare_t>
      template<class Expr>
      typename image_used_by_neighbor_list_traits<Expr>::type
      get_image_used_by_neighbor_list(Expr&& x);

      template<class Expr>
      typename image_used_by_neighbor_list_traits<Expr>::zip_image_type
      get_image_used_by_neighbor_list2(Expr&& x);

      template <class I>
      using null_image_t = constant_image<typename std::decay<I>::type::domain_type,
                                          dontcare_t>;

      /******************************************/
      /****          Implementation          ****/
      /******************************************/



      template <class Expr, int k>
      struct expression_traits< Expr, tag::image_call_p<k> >
      {
        typedef mpl::int_<k+1> number_of_images;
      };

      template <class Expr, int i>
      struct image_usage_traits<Expr, i, tag::image_call_p<i> >
      {
        typedef typename proto::result_of::value<Expr>::type type_;
        typedef typename std::remove_reference<type_>::type::type & type;

        //typedef typename std::add_lvalue_reference<typename std::remove_reference<type_>::type::type>::type type;
        //        typedef typename std::remove_reference<type_>::type::type & type;
        typedef mpl::true_         point;
        typedef mpl::false_        neighbor;
      };

      template <class Expr, int k>
      struct expression_traits< Expr, tag::image_call_n<k> >
      {
        typedef mpl::int_<k+1> number_of_images;
      };

      template <class Expr, int i>
      struct image_usage_traits<Expr, i, tag::image_call_n<i> >
      {
        typedef typename proto::result_of::value<Expr>::type type_;
        typedef typename std::remove_reference<type_>::type::type & type;

        //typedef typename std::add_lvalue_reference<typename std::remove_reference<type_>::type::type>::type type;

        typedef mpl::false_        point;
        typedef mpl::true_         neighbor;
      };

      template <class Expr>
      struct expression_traits< Expr, proto::tag::terminal>
      {
        typedef mpl::int_<0> number_of_images;
      };

      template <class Expr, int i>
      struct image_usage_traits<Expr, i, proto::tag::terminal>
      {
        typedef mpl::void_         type;
        typedef mpl::false_        point;
        typedef mpl::false_        neighbor;
      };

      template <class Expr, int i, int k>
      struct image_usage_traits<Expr, i, tag::image_call_n<k> >
      {
        typedef mpl::void_         type;
        typedef mpl::false_        point;
        typedef mpl::false_        neighbor;
      };

      template <class Expr, int i, int k>
      struct image_usage_traits<Expr, i, tag::image_call_p<k> >
      {
        typedef mpl::void_         type;
        typedef mpl::false_        point;
        typedef mpl::false_        neighbor;
      };

      template <class Expr, class... ChildTrait>
      struct expression_traits_base
      {
        typedef typename mpl::deref<
          typename mpl::max_element<
            mpl::list<typename ChildTrait::number_of_images...> >
          ::type>
        ::type number_of_images;
      };

      template <class Expr, class IntSeq>
      struct expression_traits_base_;

      template <class Expr, int... k>
      struct expression_traits_base_<Expr, intseq<k...> >
        : expression_traits_base<Expr, expression_traits<typename proto::result_of::child_c<Expr, k>::type>...>
      {
      };

      template <class Expr, class Tag>
      struct expression_traits
        : expression_traits_base_<Expr, typename int_list_seq<proto::arity_of<Expr>::value>::type >
      {
      };


      template <class Expr, int i, class IntSeq>
      struct image_usage_traits_base_;

      template <class Expr, class... ChildTrait>
      struct image_usage_traits_base
      {
        typedef typename mpl::fold<
          mpl::list<typename ChildTrait::type...>,
          mpl::void_,
          mpl::if_<boost::is_same<mpl::_1, mpl::void_>, mpl::_2, mpl::_1>
          >::type type;

        typedef typename mpl::fold<
          mpl::list<typename ChildTrait::point...>,
          mpl::false_,
          mpl::or_<mpl::_1,mpl::_2> >::type point;

        typedef typename mpl::fold<
          mpl::list<typename ChildTrait::neighbor...>,
          mpl::false_,
          mpl::or_<mpl::_1,mpl::_2> >::type neighbor;
      };


      template <class Expr, int i, int... k>
      struct image_usage_traits_base_<Expr, i, intseq<k...> >
        : image_usage_traits_base<Expr, image_usage_traits<typename proto::result_of::child_c<Expr, k>::type, i>...>
      {
      };

      template <class Expr, int i, class Tag>
      struct image_usage_traits
        : image_usage_traits_base_<Expr, i, typename int_list_seq<proto::arity_of<Expr>::value>::type >
      {
      };

      /*******    image_list  traits ***********/

      template <class Expr,
                class Seq = typename int_list_seq<number_of_images<Expr>::value>::type>
      struct image_list_traits_base;

      template <class Expr, int... k>
      struct image_list_traits_base<Expr, intseq<k...> >
      {
        typedef std::tuple<typename image_usage_traits<Expr, k>::type...>  type;
        typedef zip_image<typename image_usage_traits<Expr, k>::type...>   zip_image_type;
      };

      template <class Expr>
      struct image_list_traits : image_list_traits_base<Expr>
      {
      };

      template <class Expr,
                class Ses = typename int_list_seq<number_of_images<Expr>::value>::type>
      struct image_used_by_neighbor_list_traits_base;

      template <class Expr, int... k>
      struct image_used_by_neighbor_list_traits_base<Expr, intseq<k...> >
      {
        typedef null_image_t<typename image_usage_traits<Expr, 0>::type> null_image_type;

        typedef std::tuple<
          typename std::conditional<image_usage_traits<Expr, k>::neighbor::value,
                                    typename image_usage_traits<Expr, k>::type,
                                    null_image_type>::type ...
          > type;

        typedef zip_image<
          typename std::conditional<image_usage_traits<Expr, k>::neighbor::value,
                                    typename image_usage_traits<Expr, k>::type,
                                    null_image_type>::type ...
          > zip_image_type;
      };

      template <class Expr>
      struct image_used_by_neighbor_list_traits : image_used_by_neighbor_list_traits_base<Expr>
      {
      };


      template <int k>
      struct image_getter :
        proto::or_<
        proto::when<proto::basic_expr<tag::image_call_p<k>, proto::term<proto::_> >, proto::_value>,
        proto::when<proto::basic_expr<tag::image_call_n<k>, proto::term<proto::_> >, proto::_value>,
        proto::when<proto::basic_expr<proto::_, proto::term<proto::_> >, mpl::void_ ()>,
        proto::when<
          proto::if_<
            mpl::and_<
              mpl::greater<proto::arity_of<proto::_>, mpl::int_<0> >,
              mpl::not_<boost::is_same< image_getter<k>(proto::_child0), mpl::void_> >
              >()
            >,
          image_getter<k>(proto::_child0)
          >,
        proto::when<
          proto::if_<
            mpl::and_<
              mpl::greater<proto::arity_of<proto::_>, mpl::int_<1> >,
              mpl::not_<boost::is_same< image_getter<k>(proto::_child1), mpl::void_> >
              >()
            >,
          image_getter<k>(proto::_child1)
          >
        >
      {
      };

      template <int k>
      struct image_used_by_neighbor_getter
      {
        template <class Expr>
        typename std::enable_if<
          image_usage_traits<Expr,k>::neighbor::value,
          typename image_usage_traits<Expr,k>::type>::type
        operator()  (Expr&& x) const
        {
          return (image_getter<k> () (std::forward<Expr>(x))).get();
        }

        template <class Expr>
        typename std::enable_if<
          not image_usage_traits<Expr,k>::neighbor::value,
          null_image_t<typename image_usage_traits<Expr, 0>::type>
          >::type
        operator()  (Expr&& x) const
        {
          typedef null_image_t<typename image_usage_traits<Expr, 0>::type> R;
          return R((image_getter<0> () (std::forward<Expr>(x))).get().domain(), dontcare);
        }

      };

      template <int k, class Expr>
      typename image_usage_traits<Expr, k>::type
      get_image(Expr&& x)
      {
        return image_getter<k>() (std::forward<Expr>(x)).get(); //  proto::value(*res);
      }

      template <class Expr, int... k>
      typename image_list_traits<Expr>::type
      get_image_list_helper(Expr&& x, intseq<k...>)
      {
        return std::forward_as_tuple(image_getter<k>() (std::forward<Expr>(x)) .get()...);
      }

      template <class Expr, int... k>
      typename image_list_traits<Expr>::zip_image_type
      get_image_list_helper2(Expr&& x, intseq<k...>)
      {
        return imzip(image_getter<k>() (std::forward<Expr>(x)).get()...);
      }

      template <class Expr>
      typename image_list_traits<Expr>::type
      get_image_list(Expr&& x)
      {
        constexpr int K = expression_traits<Expr>::number_of_images::value;
        return get_image_list_helper(std::forward<Expr>(x), typename int_list_seq<K>::type ());
      }

      template <class Expr>
      typename image_list_traits<Expr>::zip_image_type
      get_image_list2(Expr&& x)
      {
        constexpr int K = expression_traits<Expr>::number_of_images::value;
        return get_image_list_helper2(std::forward<Expr>(x), typename int_list_seq<K>::type ());
      }


      template <class Expr, int... k>
      typename image_used_by_neighbor_list_traits<Expr>::type
      get_image_used_by_neighbor_list_helper(Expr&& x, intseq<k...>)
      {
        return std::forward_as_tuple(image_used_by_neighbor_getter<k>() (std::forward<Expr>(x))...);
      }

      template <class Expr>
      typename image_used_by_neighbor_list_traits<Expr>::type
      get_image_used_by_neighbor_list(Expr&& x)
      {
        constexpr int K = expression_traits<Expr>::number_of_images::value;
        return get_image_used_by_neighbor_list_helper(std::forward<Expr>(x), typename int_list_seq<K>::type ());
      }

      template <class Expr, int... k>
      typename image_used_by_neighbor_list_traits<Expr>::zip_image_type
      get_image_used_by_neighbor_list_helper2(Expr&& x, intseq<k...>)
      {
        return imzip(image_used_by_neighbor_getter<k>() (std::forward<Expr>(x))...);
      }

      template <class Expr>
      typename image_used_by_neighbor_list_traits<Expr>::zip_image_type
      get_image_used_by_neighbor_list2(Expr&& x)
      {
        constexpr int K = expression_traits<Expr>::number_of_images::value;
        return get_image_used_by_neighbor_list_helper2(std::forward<Expr>(x), typename int_list_seq<K>::type ());
      }



    } // end of namespace mln::kernel::details
  } // end of namespace mln::kernel
} // end of namespace mln

#endif //!MLN_KERNELV2_DETAILS_EXPRESSIONS_TRAITS_HPP
