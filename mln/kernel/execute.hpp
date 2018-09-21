#ifndef MLN_KERNEL_EXECUTE_HPP
# define MLN_KERNEL_EXECUTE_HPP

# include <mln/core/neighborhood/neighborhood.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/image/constant_image.hpp>
# include <mln/core/image/morphers/morpher_base.hpp>
# include <mln/core/dontcare.hpp>
# include <mln/core/internal/intseq.hpp>
# include <mln/kernel/context.hpp>
# include <mln/kernel/intro.hpp>

# include <boost/fusion/include/as_list.hpp>

namespace mln
{

  namespace kernel
  {

    ///
    /// \brief Execute a kernel expression
    ///
    ///
    template <class Expr, class N, class... I>
    void
    execute(const Expr& expr, const Neighborhood<N>& nbh, I&&... images);



    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace impl
    {

      template <class ExprList, class Context, int... K>
      void
      exec_nbh(const ExprList& elist, Context& ctx, const intseq<K...>&)
      {
        // eval with the const version of the variable enables
        // gcc to perform a better SRA since ctx does not need to
        // live in memory anymore.
        namespace bf = boost::fusion;

        const Context ctx_ = ctx;
        dontcare( (bf::at_c<K>(ctx.m_accus)
                   .take(proto::eval(proto::right(bf::at_c<K>(elist)), ctx_)), true) ... );
      }

      template <class AccuList, int... K>
      void
      init_nbh(AccuList& accus, const intseq<K...>&)
      {
        dontcare( (boost::fusion::at_c<K>(accus).init(), true) ... );
      }

      template <class Context>
      struct eval_in
      {
        template <class> struct result;

        template <class This, class Expr>
        struct result<This(Expr)> : boost::proto::functional::eval::result<This(Expr,Context)>
        {
        };


        template<typename Expr>
        typename boost::proto::result_of::eval<Expr, Context>::type
        operator()(Expr expr) const;

      };

      template <class Context>
      struct meta_make_accu
      {
        template <class> struct result;

        template <class F, class Expr>
        struct result<F(Expr)>
        {
          typedef typename proto::result_of::child_c<Expr, 0>::type	A0;
	  typedef typename proto::result_of::value<A0>::value_type	feature;
          typedef typename proto::result_of::child_c<Expr, 1>::type	child;
          typedef typename proto::result_of::eval<typename std::remove_reference<child>::type, Context>::type V;
          typedef typename feature::template apply<typename std::decay<V>::type>::type  type;
        };

        template <class Expr>
        typename result< void(Expr)>::type
        operator() (const Expr& x) const
        {
          typedef typename proto::result_of::child_c<Expr, 1>::type child;
          typedef typename proto::result_of::eval<typename std::remove_reference<child>::type, Context>::type V;

          return proto::value(x.child0).template make<typename std::decay<V>::type>() ;
        }
      };

      template <class Expr, int K, class ImageTuple>
      struct get_used_image_tuple_helper
      {
	typedef typename std::remove_reference<
	  typename std::tuple_element<0, ImageTuple>::type
	  >::type first_image;

	typedef typename first_image::domain_type domain_t;

	static constexpr bool used = get_image_usage<Expr, K>::with_n::value;

	template <bool used, class Dummy = void>
	struct foo;

	template <class Dummy>
	struct foo<true, Dummy>
	{
	  typedef typename std::tuple_element<K, ImageTuple>::type result_type;

	  result_type operator() (ImageTuple& t) const
	  {
	    return std::get<K>(t);
	  }
	};

	template <class Dummy>
	struct foo<false, Dummy>
	{
	  typedef constant_image<domain_t, dontcare_t> result_type;

	  result_type operator() (ImageTuple& t) const
	  {
	    return result_type(std::get<0>(t).domain(), dontcare);
	  }
	};

	typedef typename foo<used>::result_type result_type;

	result_type
	operator () (ImageTuple& t) const
	{
	  return foo<used> () (t);
	}
      };



      template <class Expr, int... K, class... I>
      zip_image< typename
		 get_used_image_tuple_helper<Expr, K, std::tuple<I...> >::result_type ...
		 >
      get_used_image_tuple(const Expr&, const intseq<K...>&, I&&... images)
      {
	auto t = std::forward_as_tuple(images...);

	return imzip( get_used_image_tuple_helper<Expr, K, std::tuple<I...> > () (t) ... );
      }

      template <class Pixter, class I>
      struct kernel_wrap_pixter_pixel
        : morpher_pixel_base< kernel_wrap_pixter_pixel<Pixter, I>,
                              typename Pixter::value_type>
      {
      public:
        friend struct mln::morpher_core_access;
        typedef I image_type;
        typedef typename I::reference reference;
        typedef typename I::value_type value_type;

        kernel_wrap_pixter_pixel(I* ima, const typename Pixter::reference& pix)
          : m_ima(ima), m_pix(pix)
        {
        }

        image_type& image() const { return *m_ima; }
        reference   val()   const { return m_pix.val(); }

      private:
        I*                              m_ima;
        typename Pixter::reference      m_pix;
      };


      template <class Pixter, class I>
      struct kernel_wrap_pixter :
        iterator_base< kernel_wrap_pixter<Pixter, I>,
                       kernel_wrap_pixter_pixel<Pixter, I>,
                       kernel_wrap_pixter_pixel<Pixter, I> >
      {
        kernel_wrap_pixter(Pixter* pixter, I* ima)
          : m_ima(ima), m_pixter(pixter)
        {
        }

        void init() { m_pixter->init(); }
        void next() { m_pixter->next(); }
        bool finished() const { return m_pixter->finished(); }

        kernel_wrap_pixter_pixel<Pixter, I>
        dereference() const
        {
          return kernel_wrap_pixter_pixel<Pixter, I>(m_ima, *(*m_pixter));
        }

      private:
        I*      m_ima;
        Pixter* m_pixter;
      };



      template <class Expr, class Nbh, class... I>
      void
      execute(Expr& expr, const Nbh& nbh, I&&... images)
      {


        kernel::details::transform_aggregate<> t1;
        kernel::details::retrieve_aggregate  t2;

        auto newexpr = t1(expr);
        auto subexprs_ = t2(expr, boost::fusion::nil());
        auto subexprs = boost::fusion::as_list(subexprs_);

        // Defines some helperf for aggregates
        // SubexprList: The type of the fusion list holding subexprs
        // SubExprTypeList: The type of the mpl list holding the subexprs evaluation type
        typedef decltype(subexprs) SubexprList;
        typedef meta_kernel_context<I...> meta_ctx;
        typedef eval_in<meta_ctx> Acc2Type;
        typedef typename boost::fusion::result_of::transform<const SubexprList, Acc2Type>::type SubExprTypeList_;
        typedef typename boost::fusion::result_of::as_list<SubExprTypeList_>::type SubExprTypeList;

        // accumulator list for each aggregate
        typedef typename boost::fusion::result_of::transform<const SubexprList, meta_make_accu<meta_ctx> >::type AccuList_;
        typedef typename boost::fusion::result_of::as_list<AccuList_>::type AccuList;
        AccuList_ accus_ = boost::fusion::transform(subexprs, meta_make_accu<meta_ctx> ());
        AccuList accus = boost::fusion::as_list(accus_);

        typename int_list_seq<boost::mpl::size<SubexprList>::type::value>::type iseq;

        typedef zip_image<I...> image_t;
        typedef mln_reference(image_t) V;




        // Wraps the temporary pixel iterator in a new one that
        // whose image is the tuple image where unused image on neighborhood
        // are replaced by a null
        auto z = imzip(images...);
        mln_pixter(px, z);


	auto z2 = get_used_image_tuple(expr, typename int_list_seq<sizeof... (I)>::type (), images...);
	typedef decltype(z2) image_t_2;
	typedef mln_reference(image_t_2) V2;

        kernel_wrap_pixter<decltype(px), image_t_2> tmp(&px, &z2);
        mln_iter(nx, nbh(tmp));

        // Define context
        typedef kernel::kernel_context<V, V2, SubExprTypeList, AccuList> Context;
        //typedef kernel::kernel_context<V, V, SubExprTypeList, AccuList>  Context2;

        mln_forall(px)
        {
            init_nbh(accus, iseq);

            mln_forall(nx)
            {
              Context ctx (accus, px->val(), nx->val());
              exec_nbh(subexprs, ctx, iseq);
            }

            const Context ctx (accus, px->val(), px->val());
            proto::eval(newexpr, ctx);
          }
      }

    }


    template <class Expr, class N, class... I>
    void
    execute(const Expr& expr, const Neighborhood<N>& nbh, I&&... images)
    {
      impl::execute(expr, exact(nbh), std::forward<I>(images)...);
    }
  }

} // end of namespace mln::kernel

#endif //!MLN_KERNEL_EXECUTE_HPP
