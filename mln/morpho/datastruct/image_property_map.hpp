#ifndef MLN_MORPHO_DATASTRUCT_IMAGE_PROPERTY_MAP_HPP
# define MLN_MORPHO_DATASTRUCT_IMAGE_PROPERTY_MAP_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/neighborhood/neighborhood_base.hpp>
# include <mln/core/object_wrappers.hpp>
# include <mln/morpho/datastruct/component_tree.hpp>
# include <mln/morpho/datastruct/attribute_map.hpp>
# include <mln/core/range/transform.hpp>

namespace mln
{

  namespace morpho
  {
    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map;

    // The neighborhood on the tree
    struct tree_neighb_t;

    template <class P, class AMap, class ValueMap>
    image_tree_property_map<P, AMap, ValueMap>
    make_image(const component_tree<P, AMap >& tree,
               const ValueMap& vmap);
  }


  template <class P, class AMap, class ValueMap>
  struct image_traits< morpho::image_tree_property_map<P, AMap, ValueMap> >
  {
    typedef std::true_type              accessible;
    typedef random_access_image_tag     category;
    typedef std::false_type             concrete; // depend
    typedef std::true_type              indexable;
    typedef extension::none_extension_tag         extension;
    typedef std::false_type             shallow_copy;
  };

  template <class P, class AMap, class ValueMap>
  struct image_concrete< morpho::image_tree_property_map<P,AMap,ValueMap> >
  {
    typedef property_map<morpho::component_tree<P, AMap >, typename ValueMap::value_type> NewValueMap;
    typedef morpho::image_tree_property_map<P, AMap, NewValueMap> type;
  };

  template <class P, class AMap, class ValueMap, typename V>
  struct image_ch_value< morpho::image_tree_property_map<P,AMap,ValueMap>, V>
  {
    typedef property_map<morpho::component_tree<P, AMap >, V> NewValueMap;
    typedef morpho::image_tree_property_map<P, AMap, NewValueMap> type;
  };

  namespace morpho
  {
    // Fwd
    template <class ImTreePropMap>
     struct image_tree_property_map_pixel;

    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map :
      image_base< image_tree_property_map<P, AMap, ValueMap>,
                  typename component_tree<P, AMap>::node_type,
                  typename ValueMap::value_type >
    {
    private:
      typedef component_tree<P, AMap>           tree_t;
      typedef typename tree_t::node_type        node_t;
      typedef image_tree_property_map           self_t;

      struct fun_viter;
      struct fun_cviter;
      struct fun_pixter;
      struct fun_cpixter;


    public:
      typedef typename ValueMap::value_type                     value_type;
      typedef typename ValueMap::reference                      reference;
      typedef typename ValueMap::const_reference                const_reference;
      typedef typename tree_t::node_type                        point_type;
      typedef typename tree_t::vertex_id_t                      size_type;
      typedef void                                              difference_type;
      typedef image_tree_property_map_pixel<self_t>             pixel_type;
      typedef image_tree_property_map_pixel<const self_t>       const_pixel_type;

      template <class ImTreePropMap>
      friend struct image_tree_property_map_pixel;

      template <class, class, class>
      friend struct image_tree_property_map;


      // Ranges & iterators
      typedef typename tree_t::node_range                                       domain_type;

      typedef transformed_range<domain_type, fun_viter>      value_range;
      typedef transformed_range<domain_type, fun_cviter>     const_value_range;
      typedef transformed_range<domain_type, fun_pixter>     pixel_range;
      typedef transformed_range<domain_type, fun_cpixter>    const_pixel_range;

      image_tree_property_map(const tree_t& tree, const ValueMap& vmap);

      /// \{
      /// \brief Initializer constructors
      template <class OtherVMap>
      image_tree_property_map(const image_tree_property_map<P, AMap, OtherVMap>& other, mln::init);

      template <class OtherVMap>
      image_tree_property_map(const image_tree_property_map<P, AMap, OtherVMap>& other, const value_type& v);
      /// \}


      // Acces
      reference         operator() (const node_t& node);
      const_reference   operator() (const node_t& node) const;
      reference         at (const node_t& node);
      const_reference   at (const node_t& node) const;
      reference         operator[] (size_type index);
      const_reference   operator[] (size_type index) const;
      pixel_type        pixel(const node_t& node);
      const_pixel_type  pixel(const node_t& node) const;

      node_t          point_at_index(size_type index) const;
      size_type       index_of_point(const node_t& node) const;


      /// This function does not make sense on a Tree Property Map.
      /// It's here to fulfill the concept requirement.
      difference_type delta_index(const node_t& node) const;
      void            reindex(size_type index);

      // Range & iterator
      domain_type               domain() const;
      value_range               values();
      const_value_range         values() const;
      pixel_range               pixels();
      const_pixel_range         pixels() const;

      // Specific.
      const ValueMap&           get_vmap() const;
      ValueMap&                 get_vmap();

    private:
      tree_t           m_tree;
      ValueMap         m_vmap;
    };

    /******************************/
    /*** Implementation         ***/
    /******************************/

    template <class ImTreePropMap>
    struct image_tree_property_map_pixel
      : Pixel< image_tree_property_map_pixel<ImTreePropMap> >
    {
      typedef ImTreePropMap                             image_type;
      typedef mln_value(ImTreePropMap)                  value_type;
      typedef mln_reference(ImTreePropMap)              reference;
      typedef mln_point(ImTreePropMap)                  point_type;
      typedef mln_point(ImTreePropMap)                  site_type;
      typedef typename ImTreePropMap::size_type         size_type;

      friend typename ImTreePropMap::const_pixel_type;

      image_tree_property_map_pixel() = default;
      image_tree_property_map_pixel(ImTreePropMap* ima, size_type node_id)
        : m_ima(ima), m_index(node_id)
      {
      }

      image_tree_property_map_pixel(const typename ImTreePropMap::pixel_type& other)
        : m_ima(other.m_ima), m_index(other.m_index)
      {
      }

      reference val() const
      {
        return m_ima->m_vmap[m_index];
      }

      point_type point() const
      {
        return m_ima->m_tree.get_node(m_index);
      }

      site_type site() const
      {
        return m_ima->m_tree.get_node(m_index);
      }

      image_type& image() const
      {
        return *m_ima;
      }

      size_type index() const
      {
        return m_index;
      }

    private:
      image_type*               m_ima;
      size_type                 m_index;
    };

    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map<P, AMap, ValueMap>::fun_viter
    {
      reference
      operator() (const point_type& node) const
      {
        return m_this->m_vmap[node];
      }

      self_t* m_this;
    };

    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map<P, AMap, ValueMap>::fun_cviter
    {
      fun_cviter() = default;
      fun_cviter(const self_t* self) : m_this(self) {};
      fun_cviter(const fun_viter& other) : m_this(other.m_this) {}

      const_reference
      operator() (const point_type& node) const
      {
        return m_this->m_vmap[node];
      }

      const self_t* m_this;
    };

    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map<P, AMap, ValueMap>::fun_pixter
    {
      pixel_type
      operator() (const point_type& node) const
      {
        return pixel_type(m_this, node.id());
      }

      self_t* m_this;
    };

    template <class P, class AMap, class ValueMap>
    struct image_tree_property_map<P, AMap, ValueMap>::fun_cpixter
    {
      fun_cpixter() = default;
      fun_cpixter(const self_t* self) : m_this(self) {};
      fun_cpixter(const fun_pixter& other) : m_this(other.m_this) {}

      const_pixel_type
      operator() (const point_type& node) const
      {
        return const_pixel_type(m_this, node.id());
      }

      const self_t* m_this;
    };

    template <class P, class AMap, class ValueMap>
    image_tree_property_map<P, AMap, ValueMap>::image_tree_property_map(const tree_t& tree, const ValueMap& vmap)
      : m_tree(tree),
        m_vmap(vmap)
    {
    }

    template <class P, class AMap, class ValueMap>
    template <class OtherVMap>
    image_tree_property_map<P, AMap, ValueMap>::image_tree_property_map(const image_tree_property_map<P, AMap, OtherVMap>& other,
                                                                        const value_type& v)
      : m_tree(other.m_tree),
        m_vmap(m_tree, v)
    {
    }

    template <class P, class AMap, class ValueMap>
    template <class OtherVMap>
    image_tree_property_map<P, AMap, ValueMap>::image_tree_property_map(const image_tree_property_map<P, AMap, OtherVMap>& other,
                                                                        mln::init)
      : m_tree(other.m_tree),
        m_vmap(m_tree)
    {
    }


    template <class P, class AMap, class ValueMap>
    inline
    typename image_tree_property_map<P, AMap, ValueMap>::reference
    image_tree_property_map<P, AMap, ValueMap>::operator() (const node_t& node)
    {
      return m_vmap[node];
    }


    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_reference
    image_tree_property_map<P, AMap, ValueMap>::operator() (const node_t& node) const
    {
      return m_vmap[node];
    }
    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::reference
    image_tree_property_map<P, AMap, ValueMap>::at (const node_t& node)
    {
      return m_vmap[node];
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_reference
    image_tree_property_map<P, AMap, ValueMap>::at (const node_t& node) const
    {
      return m_vmap[node];
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::reference
    image_tree_property_map<P, AMap, ValueMap>::operator[] (size_type index)
    {
      return m_vmap[index];
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_reference
    image_tree_property_map<P, AMap, ValueMap>::operator[] (size_type index) const
    {
      return m_vmap[index];
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_pixel_type
    image_tree_property_map<P, AMap, ValueMap>::pixel(const node_t& node) const
    {
      return {this, node.id()};
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::pixel_type
    image_tree_property_map<P, AMap, ValueMap>::pixel(const node_t& node)
    {
      return {this, node.id()};
    }


    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::point_type
    image_tree_property_map<P, AMap, ValueMap>::point_at_index(size_type index) const
    {
      return m_tree.get_node(index);
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::size_type
    image_tree_property_map<P, AMap, ValueMap>::index_of_point(const node_t& node) const
    {
      return node.id();
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::difference_type
    image_tree_property_map<P, AMap, ValueMap>::delta_index(const node_t& node) const
    {
      (void) node;
    }

    template <class P, class AMap, class ValueMap>
    void
    image_tree_property_map<P, AMap, ValueMap>::reindex(size_type i)
    {
      (void) i;
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::domain_type
    image_tree_property_map<P, AMap, ValueMap>::domain() const
    {
      return m_tree.nodes();
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::value_range
    image_tree_property_map<P, AMap, ValueMap>::values()
    {
      fun_viter fun = { this };
      return value_range(m_tree.nodes(), fun);
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_value_range
    image_tree_property_map<P, AMap, ValueMap>::values() const
    {
      fun_cviter fun = { this };
      return const_value_range(m_tree.nodes(), fun);
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::pixel_range
    image_tree_property_map<P, AMap, ValueMap>::pixels()
    {
      fun_pixter fun = { this };
      return pixel_range(m_tree.nodes(), fun);
    }

    template <class P, class AMap, class ValueMap>
    typename image_tree_property_map<P, AMap, ValueMap>::const_pixel_range
    image_tree_property_map<P, AMap, ValueMap>::pixels() const
    {
      fun_cpixter fun = { this };
      return const_pixel_range(m_tree.nodes(), fun);
    }

    template <class P, class AMap, class ValueMap>
    const ValueMap&
    image_tree_property_map<P, AMap, ValueMap>::get_vmap() const
    {
      return m_vmap;
    }

    template <class P, class AMap, class ValueMap>
    ValueMap&
    image_tree_property_map<P, AMap, ValueMap>::get_vmap()
    {
      return m_vmap;
    }


    template <class P, class AMap, class ValueMap>
    image_tree_property_map<P, AMap, ValueMap>
    make_image(const component_tree<P, AMap>& tree,
               const ValueMap& vmap)
    {
      return image_tree_property_map<P, AMap, ValueMap>(tree, vmap);
    }


    /***********************************/
    /** Neighborhood iterators       ***/
    /***********************************/

    namespace internal
    {

      template <class Wrapper,
                typename _is_iterator = typename mln::is_a<typename Wrapper::type, Iterator>::type>
      struct tree_nbh_iter_helper;

      template <class Wrapper>
      struct tree_nbh_iter_helper<Wrapper, std::false_type>
      {
        typedef typename Wrapper::type value_type;

      protected:
        template <class U>
        const value_type& __get(const U& x) const { return x.get(); }
      };

      template <class Wrapper>
      struct tree_nbh_iter_helper<Wrapper, std::true_type>
      {
        typedef typename Wrapper::type::value_type value_type;
      protected:

        template <class U>
        typename Wrapper::type::reference __get(const U& x) const { return *(x.get()); }
      };


      template <class NodeType>
      struct tree_nbh_piter
      : tree_nbh_iter_helper<NodeType>,
        iterator_base< tree_nbh_piter<NodeType>,
                       typename tree_nbh_iter_helper<NodeType>::value_type,
                       typename tree_nbh_iter_helper<NodeType>::value_type >
      {
        typedef typename tree_nbh_iter_helper<NodeType>::value_type node_type;
        typedef typename node_type::tree_t tree_t;
        typedef node_type value_type;

        tree_nbh_piter() = default;

        tree_nbh_piter(const NodeType& x)
        : m_x (x)
        {
        }

        void init()
        {
          m_started = (this->__get(m_x).get_parent_id() == node_type::tree_t::npos());
          m_child_iter = this->__get(m_x).children().iter();
          m_child_iter.init();
        }

        bool finished() const
        {
          return m_started and m_child_iter.finished();
        }

        void next()
        {
          if (m_started)
            m_child_iter.next();
          m_started = true;
        }

        node_type dereference() const
        {
          return m_started ? (*m_child_iter) : this->__get(m_x).parent();
        }

      private:
        NodeType                                m_x; // A node or a node iterator
        bool                                    m_started;
        typename tree_t::children_iterator      m_child_iter;
      };

      template <class NodePixelType>
      struct tree_nbh_pixter
      : tree_nbh_iter_helper<NodePixelType>,
        iterator_base< tree_nbh_pixter<NodePixelType>,
                       typename tree_nbh_iter_helper<NodePixelType>::value_type,
                       typename tree_nbh_iter_helper<NodePixelType>::value_type >
      {
        typedef typename tree_nbh_iter_helper<NodePixelType>::value_type value_type;
        typedef typename value_type::point_type node_type;
        typedef typename node_type::tree_t tree_t;

        tree_nbh_pixter() = default;
        tree_nbh_pixter(const NodePixelType& x)
        : m_x (x)
        {
        }

        void init()
        {
          auto x = this->__get(m_x).point();
          m_started = (x.get_parent_id() == node_type::tree_t::npos());
          m_child_iter = x.children().iter();
          m_child_iter.init();
        }

        bool finished() const
        {
          return m_started and m_child_iter.finished();
        }

        void next()
        {
          if (m_started)
            m_child_iter.next();
          m_started = true;
        }

        value_type dereference() const
        {
          auto x = this->__get(m_x);
          return m_started ?
            value_type(& (x.image()), m_child_iter->id())  :
            value_type(& (x.image()), x.point().get_parent_id()) ;
        }

      private:
        NodePixelType                           m_x; // A node or a node iterator
        bool                                    m_started;
        typename tree_t::children_iterator      m_child_iter;
      };

    }


    // The neighborhood on the tree
    struct tree_neighb_t : neighborhood_base<tree_neighb_t,
                                             adaptative_neighborhood_tag>
    {
    public:
      template <class NodeType>
      iterator_range< internal::tree_nbh_piter< object_wrapper<NodeType> > >
      __process_point(const NodeType& node) const
      {
        typedef iterator_range< internal::tree_nbh_piter< object_wrapper<NodeType> > > R;
        return R{ { node } };
      }

      template <class NodeType>
      iterator_range< internal::tree_nbh_piter<std::reference_wrapper<const NodeType> > >
      __bind_point(NodeType& node) const
      {
        typedef iterator_range< internal::tree_nbh_piter<std::reference_wrapper<const NodeType> > > R;
        return R{ { std::cref(node) } };
      }

      template <class NodeIterator>
      iterator_range< internal::tree_nbh_piter<std::reference_wrapper<const NodeIterator> > >
      __bind_point_iterator(const NodeIterator& node_iter) const
      {
        typedef iterator_range< internal::tree_nbh_piter<std::reference_wrapper<const NodeIterator> > > R;
        return R{ { std::cref(node_iter) } };
      }

      template <class NodePixType>
      iterator_range< internal::tree_nbh_pixter<std::reference_wrapper<const NodePixType> > >
      __bind_pixel(NodePixType& node) const
      {
        typedef iterator_range< internal::tree_nbh_pixter<std::reference_wrapper<const NodePixType> > > R;
        return R{ { std::cref(node) } };
      }

      template <class NodePixIterator>
      iterator_range< internal::tree_nbh_pixter<std::reference_wrapper<const NodePixIterator> > >
      __bind_pixel_iterator(const NodePixIterator& node_iter) const
      {
        typedef iterator_range< internal::tree_nbh_pixter<std::reference_wrapper<const NodePixIterator> > > R;
        return R{ { std::cref(node_iter) } };
      }

    };

  }

}

#endif // ! MLN_MORPHO_DATASTRUCT_IMAGE_PROPERTY_MAP_HPP
