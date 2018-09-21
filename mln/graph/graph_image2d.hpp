#ifndef GRAPH_IMAGE2D_HPP
# define GRAPH_IMAGE2D_HPP

# include <mln/core/domain/box.hpp>
# include <mln/core/range/iterator_range.hpp>
# include <mln/core/image/internal/nested_loop_iterator.hpp>
# include <mln/core/range/filter.hpp>
# include <mln/core/neighb2d.hpp>

namespace mln
{

  namespace graph
  {

    namespace internal
    {

      template <typename Vtype, typename Etype, unsigned dim>
      struct undirected_graph_ndimage_data
      {
	///
	/// \param shape number of vertices by dimension
	/// \param border number of vertices in the border
	/// \param v initialization value of the vertices
	/// \param e initialization value of the edges
	undirected_graph_ndimage_data(size_t* shape, unsigned border, Vtype v = Vtype(), Etype e = Etype());
	~undirected_graph_ndimage_data();

	template <unsigned d>
	typename std::enable_if<d == dim-1>::type
	m_construct_rec(char* ptr, const Vtype& v, const Etype& e);

	template <unsigned d>
	typename std::enable_if<d != dim-1>::type
	m_construct_rec(char* ptr, const Vtype& v, const Etype& e);

	template <unsigned d>
	typename std::enable_if<d == dim-1>::type
	m_destruct_rec(char* ptr);

	template <unsigned d>
	typename std::enable_if<d != dim-1>::type
	m_destruct_rec(char* ptr);



	static constexpr std::size_t m_alignment = std::alignment_of<Vtype>::value > std::alignment_of<Etype>::value ? std::alignment_of<Vtype>::value : std::alignment_of<Etype>::value;
	typedef typename std::aligned_storage<sizeof(Vtype), m_alignment>::type AVtype;
	typedef typename std::aligned_storage<sizeof(Etype), m_alignment>::type AEtype;

	size_t m_shape[dim];
	size_t m_estrides[dim];
	size_t m_vstrides[dim];
	size_t m_nbytes;
	char*  m_buffer;
	std::allocator<Vtype> m_valloca;
	std::allocator<Etype> m_ealloca;
      };

      // FWD
      template <typename P>
      struct graph_edge_visitor_forward;

      template <typename P, size_t N>
      inline
      std::array<P, N> arr_x2(const std::array<P, N>& arr)
      {
	std::array<P, N> out;
	for (unsigned i = 0; i < N; ++i)
	  out[i] = arr[i] * 2;
	return out;
      }


    }

    template <typename Vtype, typename Etype, typename Nbh>
    struct undirected_graph_image2d
    {
    private:
      typedef mln::internal::nested_loop_iterator<
      internal::graph_edge_visitor_forward< point<short, 2> >,
      mln::internal::no_op_visitor,
      mln::internal::no_op_visitor,
      mln::internal::point_structure<short, 2>,
      mln::internal::deref_return_point_policy> edges_iterator;


      static const int N = Nbh::static_size;

      typedef iterator_range< sliding_piter<const point2d*, std::array<point2d, N> > >  nbh_range;
      typedef decltype( std::bind(&box2d::has, (const box2d*) NULL, std::placeholders::_1) ) domain_has_p_t;


    public:
      typedef strided_box<short,2>		vertices_range;
      typedef iterator_range<edges_iterator>	edges_range;
      typedef mln::filtered_range<nbh_range, domain_has_p_t> adjacency_vertex_range;
      typedef adjacency_vertex_range adjacency_edge_range;


    //typedef ...		    adjacency_vertex_range;
    //typedef ...		    adjacency_edge_range;



      undirected_graph_image2d(const box2d& domain, const Nbh& nbh, unsigned border = 3, const Vtype& v = Vtype(), const Etype& e = Etype());

      vertices_range	vertices() const;
      edges_range	edges() const;

      point2d  source(const point2d& e) const;
      point2d  target(const point2d& e) const;

      adjacency_edge_range  edges(const point2d& v) const;
      adjacency_vertex_range adjacent_vertices(const point2d& v) const;

      const Vtype& vertex(const point2d& v) const;
      Vtype& vertex(const point2d& v);
      const Etype& edge(const point2d& e) const;
      Etype& edge(const point2d& e);

      const Vtype& vertex_at(const point2d& v) const;
      Vtype& vertex_at(const point2d& v);
      const Etype& edge_at(const point2d& e) const;
      Etype& edge_at(const point2d& e);


    protected:
      static bool _is_edge(const point2d&);
      static bool _is_vertex(const point2d&);




    private:
      std::array<point2d, N> m_enbh;
      std::array<point2d, N> m_vnbh;

    protected:
      box2d		m_domain;

    private:
      std::shared_ptr< internal::undirected_graph_ndimage_data<Vtype, Etype, 2> > m_data;


      char*		m_ptr_first;  ///< Pointer to the beginning of the buffer (outside domain if there's a border)
      char*		m_ptr_pmin;   ///< Pointer to the vertex at pmin
      char*		m_ptr_origin; ///< Pointer to the vertex at (0,0)
      size_t		m_vstrides[2]; ///< Strides between a vertex line and an edge line
      size_t		m_estrides[2]; ///< Strides between an edge line and a vertex line
      size_t		m_strides[2];  ///< Strides between a vertex line and the next vertex line (ie vstride + estride)
    };

  }

}

/*
namespace boost
{

  template <typename NV, typename EV, typename Nbh>
  struct graph_traits< mln::undirected_graph_image2d<NV, EV, Nbh> >
  {
    typedef mln::point2d vertex_descriptor;
    typedef mln::point2d edge_descriptor;
    typedef boost::undirected_tag directed_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  };

}
*/

namespace mln
{
  namespace graph
  {

  namespace internal
  {

    template <typename Vtype, typename Etype, unsigned dim>
    undirected_graph_ndimage_data<Vtype, Etype, dim>::undirected_graph_ndimage_data(size_t* shape_, unsigned border, Vtype v, Etype e)
    {
      m_vstrides[dim-1] = sizeof(AVtype);
      m_estrides[dim-1] = sizeof(AEtype);
      //std::cout << "Vtype: " << sizeof(Vtype) << "-->" << sizeof(AVtype) << std::endl;
      //std::cout << "Etype: " << sizeof(Etype) << "-->" << sizeof(AEtype) << std::endl;

      m_shape[dim-1] = shape_[dim-1] + 2*border;

      for (int i = dim-2; i >= 0; --i) {
	m_shape[i] = shape_[i] + 2 * border;
	m_vstrides[i] = m_shape[i+1] * m_vstrides[i+1] + (m_shape[i+1]-1) * m_estrides[i+1];
	m_estrides[i] = (2*m_shape[i+1] - 1) * m_estrides[i+1];
      }

      m_nbytes = m_shape[0] * m_vstrides[0] + (m_shape[0]-1) * m_estrides[0];

      // Allocate data
      m_buffer = new char[m_nbytes];

      // Construct
      m_construct_rec<0>(m_buffer, v, e);
    }

    template <typename Vtype, typename Etype, unsigned dim>
    undirected_graph_ndimage_data<Vtype, Etype, dim>::~undirected_graph_ndimage_data()
    {
      m_destruct_rec<0>(m_buffer);
    }

    template <typename Vtype, typename Etype, unsigned dim>
    template <unsigned d>
    inline
    typename std::enable_if<d != dim-1>::type
    undirected_graph_ndimage_data<Vtype, Etype, dim>::m_construct_rec(char* ptr, const Vtype& v, const Etype& e)
    {
      m_construct_rec<d+1>(ptr, v, e);
      ptr += m_vstrides[d];
      for (unsigned k = 1; k < m_shape[d]; ++k)
	{
	  {
	    AEtype* tmp = (AEtype*) ptr;
	    AEtype* end = (AEtype*) (ptr+m_estrides[d]);
	    for (;tmp != end; ++tmp)
	      m_ealloca.construct(reinterpret_cast<Etype*>(tmp), e);
	  }
	  ptr += m_estrides[d];
	  //std::cout << k << "b" << std::endl;
	  m_construct_rec<d+1>(ptr, v, e); // < vertex
	  ptr += m_vstrides[d];
	}
    }

    template <typename Vtype, typename Etype, unsigned dim>
    template <unsigned d>
    inline
    typename std::enable_if<d == dim-1>::type
    undirected_graph_ndimage_data<Vtype, Etype, dim>::m_construct_rec(char* ptr, const Vtype& v, const Etype& e)
    {
      m_valloca.construct(reinterpret_cast<Vtype*>(ptr), v);
      ptr += sizeof(AVtype);
      for (unsigned k = 1; k < m_shape[d]; ++k)
	{
	  m_ealloca.construct(reinterpret_cast<Etype*>(ptr), e);
	  ptr += sizeof(AEtype);
	  m_valloca.construct(reinterpret_cast<Vtype*>(ptr), v);
	  ptr += sizeof(AVtype);
	}
    }

    template <typename Vtype, typename Etype, unsigned dim>
    template <unsigned d>
    inline
    typename std::enable_if<d != dim-1>::type
    undirected_graph_ndimage_data<Vtype, Etype, dim>::m_destruct_rec(char* ptr)
    {
      m_destruct_rec<d+1>(ptr);
      ptr += m_vstrides[d];
      for (unsigned k = 1; k < m_shape[d]; ++k)
	{
	  for (unsigned i = 0; i < m_estrides[d]; i += sizeof(AEtype))
	    m_ealloca.destroy(reinterpret_cast<Etype*>(ptr + i));

	  ptr += m_estrides[d];
	  m_destruct_rec<d+1>(ptr); // < vertex
	  ptr += m_vstrides[d];
	}
    }

    template <typename Vtype, typename Etype, unsigned dim>
    template <unsigned d>
    inline
    typename std::enable_if<d == dim-1>::type
      undirected_graph_ndimage_data<Vtype, Etype, dim>::m_destruct_rec(char* ptr)
    {
      m_valloca.destroy(reinterpret_cast<Vtype*>(ptr));
      ptr += sizeof(AVtype);
      for (unsigned k = 1; k < m_shape[d]; ++k)
	{
	  m_ealloca.destroy(reinterpret_cast<Etype*>(ptr));
	  ptr += sizeof(AEtype);
	  m_valloca.destroy(reinterpret_cast<Vtype*>(ptr));
	  ptr += sizeof(AVtype);
	}
    }


    template <typename P>
    struct graph_edge_visitor_forward
    {
      typedef P point_type;
      enum { dim = P::ndim };

      graph_edge_visitor_forward() : pmin_ (), pmax_ () {}
      graph_edge_visitor_forward(const P& pmin, const P& pmax) :
	pmin_ (pmin), pmax_ (pmax) {}

      void  initialize(P& point) const
      {
	point = pmin_;
	point[dim-1] += 1;
      }

      template <size_t n> void  init(P& point) const
      {
	if (n == 0 or point[n-1] % 2 != 0) // on a edge
	  point[n] = pmin_[n];
	else
	  point[n] = pmin_[n] + 1;
      }

      template <size_t n> void  next(P& point) const
      {
	if (n == 0 or point[n-1] % 2 != 0)
	  point[n] += 1;
	else
	  point[n] += 2;
      }
      template <size_t n> bool  finished(const P& point) const
      {
	return point[n] >= pmax_[n];
      }

    private:
      P pmin_;
      P pmax_;
    };


  }

  template <typename Vtype, typename Etype, typename Nbh>
  undirected_graph_image2d<Vtype, Etype, Nbh>::undirected_graph_image2d(const box2d& domain, const Nbh& nbh, unsigned border, const Vtype& v, const Etype& e)
  {
    size_t shp[2];
    shp[0] = (domain.pmax[0] - domain.pmin[0]);
    shp[1] = (domain.pmax[1] - domain.pmin[1]);

    m_domain.pmin = domain.pmin * 2;
    m_domain.pmax = domain.pmax * 2 - 1;

    m_data.reset(new internal::undirected_graph_ndimage_data<Vtype, Etype, 2>(&shp[0], border, v, e));
    m_vstrides[0] = m_data->m_vstrides[0];
    m_vstrides[1] = m_data->m_vstrides[1];
    m_estrides[0] = m_data->m_estrides[0];
    m_estrides[1] = m_data->m_estrides[1];
    m_strides[0] = m_data->m_estrides[0] + m_data->m_vstrides[0];
    m_strides[1] = m_data->m_estrides[1] + m_data->m_vstrides[1];


    m_ptr_first = m_data->m_buffer;
    m_ptr_pmin = m_ptr_first;
    m_ptr_pmin += border * m_vstrides[0] + (border - 1) * m_estrides[0];
    m_ptr_pmin += border * m_vstrides[1] + (border - 1) * m_estrides[1];
    m_ptr_origin = m_ptr_pmin;
    m_ptr_origin -= (domain.pmin[0]) * m_strides[0];
    m_ptr_origin -= (domain.pmin[1]) * m_strides[1];

    m_enbh = nbh.dpoints;
    m_vnbh = internal::arr_x2(nbh.dpoints);
  }


  template <typename Vtype, typename Etype, typename Nbh>
  inline
  point2d
  undirected_graph_image2d<Vtype, Etype, Nbh>::source(const point2d& e) const
  {
    mln_precondition(m_domain.has(e));
    mln_precondition(_is_edge(e));
    int x = e[0], y = e[1];
    typedef point2d::value_type T;

    if (x % 2 == 1)
      if (y % 2 == 1) // diagonal edge
	return point2d{(T) (x-1),(T) (y-1)};
      else
	return point2d{(T) (x-1),(T) y};
    else
      return point2d{(T) x,(T) (y-1)};
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  point2d
  undirected_graph_image2d<Vtype, Etype, Nbh>::target(const point2d& e) const
  {
    mln_precondition(m_domain.has(e));
    mln_precondition(_is_edge(e));
    int x = e[0], y = e[1];

    typedef point2d::value_type T;

    if (x % 2 == 1)
      if (y % 2 == 1) // diagonal edge
	return point2d{(T) (x+1),(T) (y+1)};
      else
	return point2d{(T) (x+1), (T) y};
    else
      return point2d{(T) x,(T) (y+1)};
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  bool
  undirected_graph_image2d<Vtype, Etype, Nbh>::_is_vertex(const point2d& v)
  {
    return v[0] % 2 == 0 and v[1] % 2 == 0;
  }




  template <typename Vtype, typename Etype, typename Nbh>
  inline
  bool
  undirected_graph_image2d<Vtype, Etype, Nbh>::_is_edge(const point2d& e)
  {
    return !_is_vertex(e);
  }


  template <typename Vtype, typename Etype, typename Nbh>
  inline
  Vtype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::vertex_at(const point2d& v)
  {
    mln_precondition(_is_vertex(v));
    char* ptr = m_ptr_origin;
    ptr += (v[0] / 2) * m_strides[0] + (v[1] / 2) * m_strides[1];
    return * (Vtype*) ptr;
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  const Vtype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::vertex_at(const point2d& v) const
  {
    mln_precondition(_is_vertex(v));
    const char* ptr = m_ptr_origin;
    ptr += (v[0] / 2) * m_strides[0] + (v[1]/2) * m_strides[1];
    return * (const Vtype*) ptr;
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  Etype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::edge_at(const point2d& e)
  {
    mln_precondition(_is_edge(e));
    char* ptr = m_ptr_origin;

    size_t s0 = (e[0] < 0) ? m_estrides[0] : m_vstrides[0];
    ptr += (e[0] / 2) * m_strides[0] + (e[0] % 2) * s0;

    if (e[0] % 2 == 0) {
      size_t s1 = (e[1] < 0) ? m_estrides[1] : m_vstrides[1];
      ptr += (e[1] / 2) * m_strides[1] + (e[1] % 2) * s1;
    } else {
      ptr += e[1] * m_estrides[1];
    }

    return * (Etype*) ptr;
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  const Etype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::edge_at(const point2d& e) const
  {
    mln_precondition(_is_edge(e));
    const char* ptr = m_ptr_origin;


    size_t s0 = (e[0] < 0) ? m_estrides[0] : m_vstrides[0];
    ptr += (e[0] / 2) * m_strides[0] + (e[0] % 2) * s0;

    if (e[0] % 2 == 0) {
      size_t s1 = (e[1] < 0) ? m_estrides[1] : m_vstrides[1];
      ptr += (e[1] / 2) * m_strides[1] + (e[1] % 2) * s1;
    } else {
      ptr += e[1] * m_estrides[1];
    }

    return * (const Etype*) ptr;
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  Vtype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::vertex(const point2d& v)
  {
    mln_precondition(m_domain.has(v));
    return vertex_at(v);
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  const Vtype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::vertex(const point2d& v) const
  {
    mln_precondition(m_domain.has(v));
    return vertex_at(v);
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  Etype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::edge(const point2d& e)
  {
    mln_precondition(m_domain.has(e));
    return edge_at(e);
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  const Etype&
  undirected_graph_image2d<Vtype, Etype, Nbh>::edge(const point2d& e) const
  {
    mln_precondition(m_domain.has(e));
    return edge_at(e);
  }

  template <typename Vtype, typename Etype, typename Nbh>
  inline
  typename undirected_graph_image2d<Vtype, Etype, Nbh>::vertices_range
  undirected_graph_image2d<Vtype, Etype, Nbh>::vertices() const
  {
    return strided_box<short, 2>(m_domain.pmin, m_domain.pmax, point2d{2,2});
  }


  template <typename Vtype, typename Etype, typename Nbh>
  inline
  typename undirected_graph_image2d<Vtype, Etype, Nbh>::edges_range
  undirected_graph_image2d<Vtype, Etype, Nbh>::edges() const
  {
    return make_iterator_range( edges_iterator( mln::internal::point_structure<short, 2> (),
						internal::graph_edge_visitor_forward<point2d> (m_domain.pmin, m_domain.pmax),
						mln::internal::no_op_visitor (),
						mln::internal::no_op_visitor () ));
  }

    template <typename Vtype, typename Etype, typename Nbh>
    inline
    typename undirected_graph_image2d<Vtype, Etype, Nbh>::adjacency_edge_range
    undirected_graph_image2d<Vtype, Etype, Nbh>::edges(const point2d& v) const
    {
      mln_precondition(_is_vertex(v));
      sliding_piter<const point2d*, std::array<point2d, N> > eit(&v, m_enbh);
      auto myrng = make_iterator_range(eit);
      return rng::filter(std::move(myrng), std::bind(&box2d::has, &m_domain, std::placeholders::_1));
    }

    template <typename Vtype, typename Etype, typename Nbh>
    inline
    typename undirected_graph_image2d<Vtype, Etype, Nbh>::adjacency_edge_range
    undirected_graph_image2d<Vtype, Etype, Nbh>::adjacent_vertices(const point2d& v) const
    {
      mln_precondition(_is_vertex(v));
      sliding_piter<const point2d*, std::array<point2d, N> > vit(&v, m_vnbh);
      auto myrng = make_iterator_range(vit);
      return rng::filter(std::move(myrng), std::bind(&box2d::has, &m_domain, std::placeholders::_1));
    }


  }

}

#endif // ! GRAPH_IMAGE2D_HPP
