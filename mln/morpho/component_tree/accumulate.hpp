#ifndef MLN_MORPHO_COMPONENT_TREE_ACCUMULATE_HPP
# define MLN_MORPHO_COMPONENT_TREE_ACCUMULATE_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/trace.hpp>

# include <mln/accu/accumulator.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/datastruct/attribute_map.hpp>

namespace mln
{

  namespace morpho
  {

    template <class P, class AssociativeMap, class I>
    property_map< component_tree<P, AssociativeMap>, mln_value(I) >
    make_attribute_map_from_image(const component_tree<P, AssociativeMap>& ctree,
				  const Image<I>& ima);

    template <class P, class AssociativeMap, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, P>::type >
    accumulate(const component_tree<P, AssociativeMap>& ctree,
	       const AccumulatorLike<Accu>& accu,
	       bool require_ordering = false);

    template <class P, class AssociativeMap, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, P>::type >
    accumulate_proper(const component_tree<P, AssociativeMap>& ctree,
                      const AccumulatorLike<Accu>& accu);


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_value(I)>::type >
    vaccumulate(const component_tree<P, AssociativeMap>& ctree,
                const Image<I>& ima,
                const AccumulatorLike<Accu>& accu,
                bool require_ordering = false);

    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
                  typename accu::result_of<Accu, mln_value(I)>::type >
    vaccumulate_proper(const component_tree<P, AssociativeMap>& ctree,
                       const Image<I>& ima,
                       const AccumulatorLike<Accu>& accu);


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_point(I)>::type >
    paccumulate(const component_tree<P, AssociativeMap>& ctree,
                const Image<I>& ima,
                const AccumulatorLike<Accu>& accu,
                bool require_ordering = false);




    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_point(I)>::type >
    paccumulate_proper(const component_tree<P, AssociativeMap>& ctree,
		       const Image<I>& ima,
		       const AccumulatorLike<Accu>& accu);


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_cpixel(I)>::type >
    pixaccumulate(const component_tree<P, AssociativeMap>& ctree,
		  const Image<I>& ima,
		  const AccumulatorLike<Accu>& accu,
		  bool require_ordering = true);

    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_cpixel(I)>::type >
    pixaccumulate_proper(const component_tree<P, AssociativeMap>& ctree,
			 const Image<I>& ima,
			 const AccumulatorLike<Accu>& accu);



    /********************************/
    /*** Implementation		  ***/
    /********************************/

    namespace impl
    {

      template <class Tree, class I, class Accu, class Pmap, class PixFunctor>
      void
      accumulate_proper_index(const Tree& tree, const I& ima, Accu acc, Pmap& vmap,
			      PixFunctor pixfun)
      {
	acc.init();

	// 1. Accumulate per-node
	property_map<Tree, Accu> accmap(tree, acc);
	mln_pixter(px, ima);
	mln_forall(px)
	  accmap[ tree.get_node_id(px->index()) ].take( pixfun(*px) );

	// 2. set value
	mln_foreach(auto x, tree.nodes())
	  vmap[x.id()] = accmap[x.id()].to_result();
      }


      template <class Tree, class I, class Accu, class Pmap, class PixFunctor>
      void
      accumulate_index(const Tree& tree, const I& ima, Accu acc, Pmap& vmap,
		       PixFunctor pixfun)
      {
	acc.init();
	property_map<Tree, Accu> accmap(tree, acc);

	// 1. Accumulate per-node
	mln_pixter(px, ima);
	mln_forall(px)
	  accmap[ tree.get_node_id(px->index()) ].take( pixfun(*px) );

	// 2. Transmit to parent & set value
	mln_reverse_foreach(auto x, tree.nodes())
	  {
	    vmap[x.id()] = accmap[x.id()].to_result();
	    accmap[x.get_parent_id()].take(accmap[x.id()]);
	  }
        vmap[Tree::npos()] = accmap[Tree::npos()].to_result();
      }

      template <class Tree, class I, class Accu, class Pmap, class IFunctor>
      void
      accumulate_ordered_index(const Tree& tree, const I& ima, Accu acc, Pmap& vmap,
			       IFunctor ifun)
      {
	(void) ima;
	acc.init();
	property_map<Tree, Accu> accmap(tree, acc);

	// Accumulate and transmit and set
	mln_reverse_foreach(auto x, tree.pset())
	  {
	    auto node = tree.get_node_at(x);
	    accmap[ node.id() ].take( ifun(x) );

	    // Transmit & set if last
            // WARning: Does not work if a node is empty
	    // if (node.first_point() == x)
	    //   {
	    //     vmap[node.id()] = accmap[ node.id() ].to_result();
	    //     accmap[ node.get_parent_id() ].take(accmap[node.id()]);
	    //   }
	  }

        // 2. Transmit to parent & set value
	mln_reverse_foreach(auto x, tree.nodes())
	  {
	    vmap[x.id()] = accmap[x.id()].to_result();
	    accmap[x.get_parent_id()].take(accmap[x.id()]);
	  }
        vmap[Tree::npos()] = accmap[Tree::npos()].to_result();
      }
    }

    template <class P, class AssociativeMap,  class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, P>::type >
    accumulate(const component_tree<P, AssociativeMap>& tree,
		const AccumulatorLike<Accu>& accu_,
		bool require_order)
    {
      mln_entering("mln::morpho::accumulate");
      auto acc = accu::make_accumulator(exact(accu_), P ());
      typedef typename accu::result_of<Accu, P>::type R;
      typedef component_tree<P, AssociativeMap> Tree;

      property_map<Tree, R> vmap(tree);
      if (not require_order) {
        auto pixfunctor = [](const mln_cpixel(AssociativeMap)& px) { return px.index(); };
        impl::accumulate_index(tree, tree._get_data()->m_pmap, acc, vmap, pixfunctor);
      } else {
        auto ifunctor = [](P i) { return i; };
        impl::accumulate_ordered_index(tree, tree._get_data()->m_pmap, acc, vmap, ifunctor);
      }
      mln_exiting();
      return vmap;
    }

    template <class P, class AssociativeMap,  class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, P>::type >
    accumulate_proper(const component_tree<P, AssociativeMap>& tree,
                      const AccumulatorLike<Accu>& accu_)
    {
      mln_entering("mln::morpho::accumulate");
      auto acc = accu::make_accumulator(exact(accu_), P ());
      typedef typename accu::result_of<Accu, P>::type R;
      typedef component_tree<P, AssociativeMap> Tree;

      property_map<Tree, R> vmap(tree);

      auto pixfunctor = [](const mln_cpixel(AssociativeMap)& px) { return px.index(); };
      impl::accumulate_proper_index(tree, tree._get_data()->m_pmap, acc, vmap, pixfunctor);
      mln_exiting();
      return vmap;
    }


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_value(I)>::type >
    vaccumulate(const component_tree<P, AssociativeMap>& tree,
		const Image<I>& ima_,
		const AccumulatorLike<Accu>& accu_,
		bool require_order)
    {
      mln_entering("mln::morpho::vaccumulate");
      auto acc = accu::make_accumulator(exact(accu_), mln_value(I) ());
      typedef typename accu::result_of<Accu, mln_value(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;

      property_map<Tree, R> vmap(tree);
      const I& ima = exact(ima_);

      if (not require_order) {
	auto pixfunctor = [](const mln_cpixel(I)& px) { return px.val(); };
	impl::accumulate_index(tree, exact(ima), acc, vmap, pixfunctor);
      } else {
	auto ifunctor = [&ima](typename I::size_type i) { return ima[i]; };
	impl::accumulate_ordered_index(tree, exact(ima), acc, vmap, ifunctor);
      }

      mln_exiting();
      return vmap;
    }

    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_value(I)>::type >
    vaccumulate_proper(const component_tree<P, AssociativeMap>& tree,
		       const Image<I>& ima,
		       const AccumulatorLike<Accu>& accu_)
    {
      mln_entering("mln::morpho::vaccumulate_proper");
      auto acc = accu::make_accumulator(exact(accu_), mln_value(I) ());
      typedef typename accu::result_of<Accu, mln_value(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;

      property_map<Tree, R> vmap(tree);

      auto pixfunctor = [](const mln_cpixel(I)& px) { return px.val(); };
      impl::accumulate_proper_index(tree, exact(ima), acc, vmap, pixfunctor);

      mln_exiting();
      return vmap;
    }


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_point(I)>::type >
    paccumulate(const component_tree<P, AssociativeMap>& tree,
		const Image<I>& ima_,
		const AccumulatorLike<Accu>& accu_,
		bool require_order)
    {
      mln_entering("mln::morpho::paccumulate");
      auto acc = accu::make_accumulator(exact(accu_), mln_point(I) ());

      typedef typename accu::result_of<Accu, mln_point(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;
      property_map<Tree, R> vmap(tree);

      const I& ima = exact(ima_);

      if (not require_order) {
	auto pixfunctor = [](const mln_cpixel(I)& px) { return px.point(); };
	impl::accumulate_index(tree, ima, acc, vmap, pixfunctor);
      } else {
	auto ifunctor = [&ima](typename I::size_type i) { return ima.point_at_index(i); };
	impl::accumulate_ordered_index(tree, ima, acc, vmap, ifunctor);
      }

      mln_exiting();
      return vmap;
    }


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_point(I)>::type >
    paccumulate_proper(const component_tree<P, AssociativeMap>& tree,
		       const Image<I>& ima,
		       const AccumulatorLike<Accu>& accu_)
    {
      mln_entering("mln::morpho::paccumulate_proper");
      auto acc = accu::make_accumulator(exact(accu_), mln_point(I) ());

      typedef typename accu::result_of<Accu, mln_point(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;
      property_map<Tree, R> vmap(tree);

      auto pixfunctor = [](const mln_cpixel(I)& px) { return px.point(); };
      impl::accumulate_proper_index(tree, exact(ima), acc, vmap, pixfunctor);

      mln_exiting();
      return vmap;
    }


    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_cpixel(I)>::type >
    pixaccumulate(const component_tree<P, AssociativeMap>& tree,
		  const Image<I>& ima_,
		  const AccumulatorLike<Accu>& accu_,
		  bool require_order)
    {
      mln_entering("mln::morpho::pixaccumulate");
      auto acc = accu::make_accumulator(exact(accu_), mln_cpixel(I) ());

      typedef typename accu::result_of<Accu, mln_pixel(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;
      property_map<Tree, R> vmap(tree);

      const I& ima = exact(ima_);

      if (not require_order) {
	auto pixfunctor = [](const mln_cpixel(I)& px) { return px; };
	impl::accumulate_index(tree, ima, acc, vmap, pixfunctor);
      } else {
	auto ifunctor = [&ima](typename I::size_type i) { return ima.pixel_at(ima.index_of_point(i)); };
	impl::accumulate_ordered_index(tree, ima, acc, vmap, ifunctor);
      }

      mln_exiting();
      return vmap;
    }

    template <class P, class AssociativeMap, class I, class Accu>
    property_map< component_tree<P, AssociativeMap>,
		  typename accu::result_of<Accu, mln_cpixel(I)>::type >
    pixaccumulate_proper(const component_tree<P, AssociativeMap>& tree,
			 const Image<I>& ima,
			 const AccumulatorLike<Accu>& accu_)
    {
      mln_entering("mln::morpho::pixaccumulate_proper");
      auto acc = accu::make_accumulator(exact(accu_), mln_cpixel(I) ());

      typedef typename accu::result_of<Accu, mln_cpixel(I)>::type R;
      typedef component_tree<P, AssociativeMap> Tree;
      property_map<Tree, R> vmap(tree);

      auto pixfunctor = [](const mln_cpixel(I)& px) { return px; };
      impl::accumulate_proper_index(tree, exact(ima), acc, vmap, pixfunctor);

      mln_exiting();
      return vmap;
    }


    template <class P, class AssociativeMap, class I>
    property_map< component_tree<P, AssociativeMap>, mln_value(I) >
    make_attribute_map_from_image(const component_tree<P, AssociativeMap>& tree,
				  const Image<I>& ima_)
    {
      mln_entering("mln::morpho::make_attribute_map_from_image");

      typedef component_tree<P, AssociativeMap> Tree;
      const I& ima = exact(ima_);

      property_map<Tree, mln_value(I)> vmap(tree);
      mln_foreach(auto x, tree.nodes())
	vmap[x.id()] = ima[  x.first_point() ];

      mln_exiting();
      return vmap;
    }

  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_ACCUMULATE_HPP
