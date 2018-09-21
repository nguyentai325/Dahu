#ifndef MLN_MORPHO_COMPONENT_TREE_IO_HPP
# define MLN_MORPHO_COMPONENT_TREE_IO_HPP

# include <mln/morpho/component_tree/component_tree.hpp>
# include <iosfwd>
# include <mln/io/imsave.hpp>
# include <mln/io/imread.hpp>

namespace mln
{

  namespace morpho
  {

    template <class P, class AssociativeMap>
    void
    save(const component_tree<P, AssociativeMap>& tree, std::ostream& os);

    template <class P, class AssociativeMap>
    void
    load(std::istream& is, component_tree<P, AssociativeMap>& tree);

    template <class P, class AssociativeMap>
    void
    load(const std::string& s, component_tree<P, AssociativeMap>& tree);


    /*************************************/
    /*****  Implementation           *****/
    /*************************************/

    /*
    // Do not extra information (version) of component tree node.
    // This is redoundant with the tree class.
    namespace boost {
      namespace serialization {

        template <class P, class AssociativeMap>
        struct implementation_level<const internal::component_tree_node<P, AssociativeMap> >
        {
          typedef mpl::integral_c_tag tag;
          typedef mpl::int_<boost::serialization::object_serializable> type;
          BOOST_STATIC_CONSTANT(int, value = type::value);
        };

        template <class P, class AssociativeMap>
        struct tracking_level<const internal::component_tree_node<P, AssociativeMap> >
        {
          typedef mpl::integral_c_tag tag;
          typedef mpl::int_<track_never> type;
          BOOST_STATIC_CONSTANT(int, value = type::value);
        };

        template <class P, class AssociativeMap>
        struct is_bitwise_serializable<const internal::component_tree_node<P, AssociativeMap> >
          : std::true_type
        {
        };

      }
    }
    */


    template <class P, class AssociativeMap>
    void
    save(const component_tree<P, AssociativeMap>& tree, std::ostream& s)
    {
      const internal::component_tree_data<P, AssociativeMap>* data = tree._get_data();

      s <<  (unsigned) data->m_nodes.size() << std::endl;
      s <<  (unsigned) data->m_S.size() << std::endl;
      s <<  (int)data->m_pset_ordered << std::endl;
      s <<  (unsigned) tree.get_root_id() << std::endl;

      //write node array
      s.write((const char*) &(data->m_nodes[0]), data->m_nodes.size() * sizeof(internal::component_tree_node));

      // write S array
      s.write((const char*) &(data->m_S[0]), data->m_S.size() * sizeof(P));

      // write pmap
      io::imsave(data->m_pmap, s);
      s.flush();
    }

    template <class P, class AssociativeMap>
    void
    load(std::istream& s, component_tree<P, AssociativeMap>& tree)
    {
      internal::component_tree_data<P, AssociativeMap>* data = tree._get_data();

      unsigned nnodes;
      unsigned npoints;
      int ordered;
      unsigned root;

      s >> nnodes;
      s >> npoints;
      s >> ordered;
      s >> root;
      s.ignore(1);

      //std::cout << nnodes << " " << npoints << " " << ordered << std::endl;
      data->m_nodes.resize(nnodes);
      data->m_S.resize(npoints);
      data->m_pset_ordered = ordered;

      // read node array
      s.read((char*) (&data->m_nodes[0]), nnodes * sizeof(internal::component_tree_node));

      // read S array
      s.read((char*) (&data->m_S[0]), npoints * sizeof(P));

      // // read pmap
      io::imread(s, data->m_pmap);

      // set the rooted
      tree = tree.get_subtree(root);
    }

    template <class P, class AssociativeMap>
    void
    load(const std::string& s, component_tree<P, AssociativeMap>& tree)
    {
      std::ifstream f(s);
      load(f, tree);
    }

    template <class P, class AssociativeMap>
    void
    save(component_tree<P, AssociativeMap>& tree, const std::string& s)
    {
      std::ofstream f(s);
      save(tree, f);
    }


  }

}

#endif // ! MLN_MORPHO_COMPONENT_TREE_IO_HPP
