#ifndef GRAPHCUT_HH
# define GRAPHCUT_HH

# include <mln/core/value/value_traits.hpp>
# include <mln/core/image/image2d.hpp>
# include <mln/graph/graph_image2d.hpp>
# include <queue>
# include <mln/core/trace.hpp>

namespace mln
{

  namespace graphcut
  {


    namespace internal
    {
      template <class DataFidelityFunction, class RegularityFunction, typename V>
      struct default_cmp;
    }


    template <typename V, class N, class DataFidelityFunction, class RegularityFunction, class Compare = internal::default_cmp<DataFidelityFunction, RegularityFunction, V> >
    double
    graphcut(const image2d<V>& ima, image2d<bool>& f, const N& nbh, DataFidelityFunction d, RegularityFunction v, const Compare& cmp = Compare ());



    /*********************/
    /** Implementation  **/
    /*********************/


    namespace internal
    {

      template <typename Vtype, typename Etype, typename Nbh>
      struct graphcut_graph_t : graph::undirected_graph_image2d<Vtype, Etype, Nbh>
      {
      private:
	typedef  graph::undirected_graph_image2d<Vtype, Etype, Nbh> base;

      public:
	enum edge_t { SOURCE_EDGE, NORMAL_EDGE, SINK_EDGE };

	struct edge_type {
	  edge_t	type;
	  point2d	p;

	  bool
	  operator < (const edge_t& other) const
	  {
	    return this->type < other.type or
	      this->type == other.type and this->p == other.p;
	  }
	};


	struct vertices_iterator : iterator_base<vertices_iterator, point2d, const point2d&>
	{
	  vertices_iterator(const sbox2d& domain, const point2d& psrc, const point2d& psink)
	    : m_src (psrc),
	      m_sink (psink),
	      m_domain_it (domain.iter())
	  {
	  }

	  void init()
	  {
	    m_status = 0;
	    m_domain_it.init();
	  }

	  void next()
	  {
	    mln_precondition(not finished());
	    if (m_status != 1) {
	      ++m_status;
	    } else {
	      m_domain_it.next();
	      if (m_domain_it.finished())
		++m_status;
	    }
	  }

	  bool finished() const
	  {
	    return m_status > 2;
	  }

	  const point2d& dereference() const
	  {
	    mln_precondition(not finished());
	    switch (m_status) {
	      case 0:  return m_src;
	      case 1:  return *m_domain_it;
	      default: return m_sink;
	    }
	  }

	private:
	  const point2d			m_src;
	  const point2d			m_sink;
	  sbox2d::iterator		m_domain_it;
	  int		m_status; // 0-> SRC, 1-> domain, 2-> sink
	};

	struct edges_iterator : iterator_base<edges_iterator, edge_type, const edge_type&>
	{
	  typedef typename graph::undirected_graph_image2d<Vtype, Etype, Nbh>::edges_range graph_edges_range;
	  typedef typename graph_edges_range::iterator graph_edges_iterator;

	  edges_iterator(const graph_edges_range& edges, const sbox2d& domain)
	    : m_edges (edges),
	      m_domain (domain),
	      m_edges_it (m_edges.iter()),
	      m_domain_it (m_domain.iter())
	  {
	  }

	  void init()
	  {
	    m_domain_it.init();
	    m_e.type = SOURCE_EDGE;
	    m_e.p = *m_domain_it;
	  }

	  void next()
	  {
	    mln_precondition(not finished());
	    switch (m_e.type)
	      {
		case SOURCE_EDGE:
		  m_domain_it.next();
		  if (!m_domain_it.finished()) {
		      m_e.p = *m_domain_it;
		  } else {
		    m_edges_it.init();
		    m_e.type = NORMAL_EDGE;
		    m_e.p = *m_edges_it;
		  }
		  break;
		case NORMAL_EDGE:
		  m_edges_it.next();
		  if (!m_edges_it.finished()) {
		    m_e.p = *m_edges_it;
		  } else {
		    m_domain_it.init();
		    m_e.type = SINK_EDGE;
		    m_e.p = *m_domain_it;
		  }
		  break;
		case SINK_EDGE:
		  m_domain_it.next();
		  m_e.p = *m_domain_it;
	      }
	  }

	  bool finished() const
	  {
	    return m_e.type == SINK_EDGE and m_domain_it.finished();
	  }

	  const edge_type& dereference() const
	  {
	    mln_precondition(not finished());
	    return m_e;
	  }

	private:
	  const graph_edges_range	m_edges;
	  const sbox2d			m_domain;
	  graph_edges_iterator		m_edges_it;
	  sbox2d::iterator		m_domain_it;
	  edge_type			m_e;
	};

	struct adjacency_vertex_iterator
	  : iterator_base<adjacency_vertex_iterator, point2d, point2d>
	{
	  typedef typename range_iterator<typename base::adjacency_vertex_range>::type inner_iterator;

	  adjacency_vertex_iterator(const inner_iterator& inner, const sbox2d& domain,
				    const point2d& psrc, const point2d& psink, const point2d& v)
	    : m_src (psrc),
	      m_sink (psink),
	      m_inner_it (inner),
	      m_domain_it (domain.iter()),
	      m_point (&v)
	  {
	  }

	  void init()
	  {
	    m_status = (int) (*m_point != m_src and *m_point != m_sink);
	    if (m_status == 0) {
	      m_domain_it.init();
	    }
	  }

	  void next()
	  {
	    mln_precondition(!finished());
	    switch (m_status)
	      {
		case 0:	m_domain_it.next();
		  if (m_domain_it.finished())
		    m_status = 4;
		  break;
		case 1: m_inner_it.init(); ++m_status; break;
		case 2:
		  m_inner_it.next();
		  if (m_inner_it.finished())
		    ++m_status;
		  break;
		default: ++m_status;
	      }
	  }

	  point2d dereference() const
	  {
	    mln_precondition(!finished());
	    switch (m_status)
	      {
		case 0:	return *m_domain_it;
		case 1: return m_src;
		case 2: return *m_inner_it;
		default: return m_sink;
	      }
	  }

	  bool finished() const
	  {
	    return m_status > 3;
	  }


	private:
	  const point2d		m_src;
	  const point2d		m_sink;
	  inner_iterator	m_inner_it;
	  sbox2d::iterator	m_domain_it;

	  const point2d*	m_point;
	  int			m_status; // 0 (SINK or SOURCE), 1 -> (INNER) SINK, 2-> (INNER) Nbh, 3 -> (INNER) SOURCE
	};


	struct adjacency_edge_iterator : iterator_base<adjacency_edge_iterator, edge_type, const edge_type&>
	{
	  typedef typename range_iterator<typename base::adjacency_edge_range>::type inner_iterator;

	  adjacency_edge_iterator(const inner_iterator& inner, const sbox2d& domain,
				  const point2d& psrc, const point2d& psink, const point2d& v)
	    : m_src (psrc),
	      m_sink (psink),
	      m_inner_it (inner),
	      m_domain_it (domain.iter()),
	      m_point (&v)
	  {
	  }

	  void init()
	  {
	    m_status = (int) (*m_point != m_src and *m_point != m_sink);
	    if (m_status == 0) {
	      m_edge.type = (*m_point == m_src) ? SOURCE_EDGE : SINK_EDGE;
	      m_domain_it.init();
	      m_edge.p = *m_domain_it;
	    } else {
	      m_edge.type = SOURCE_EDGE;
	      m_edge.p = *m_point;
	    }
	  }

	  void next()
	  {
	    mln_precondition(!finished());
	    switch (m_status)
	      {
		case 0:
		  m_domain_it.next();
		  if (m_domain_it.finished())
		    m_status = 4;
		  else
		    m_edge.p = *m_domain_it;
		  break;
		case 1:
		  m_inner_it.init();
		  ++m_status;
		  m_edge.type = NORMAL_EDGE;
		  m_edge.p = *m_inner_it;
		  break;
		case 2:
		  m_inner_it.next();
		  if (m_inner_it.finished())
		    {
		      ++m_status;
		      m_edge.type = SINK_EDGE;
		      m_edge.p = *m_point;
		    }
		  else
		    m_edge.p = *m_inner_it;
		  break;
		default: ++m_status;
	      }
	  }

	  const edge_type& dereference() const
	  {
	    mln_precondition(!finished());
	    return m_edge;
	  }

	  bool finished() const
	  {
	    return m_status > 3;
	  }


	private:
	  const point2d		m_src;
	  const point2d		m_sink;
	  inner_iterator	m_inner_it;
	  sbox2d::iterator	m_domain_it;

	  const point2d*	m_point;
	  int			m_status; // 0 (SINK or SOURCE), 1 -> (INNER) SINK, 2-> (INNER) Nbh, 3 -> (INNER) SOURCE
	  edge_type		m_edge;
	};




	typedef iterator_range<vertices_iterator> vertices_range;
	typedef iterator_range<edges_iterator>	  edges_range;
	typedef iterator_range<adjacency_vertex_iterator> adjacency_vertex_range;
	typedef iterator_range<adjacency_edge_iterator> adjacency_edge_range;

	typedef point2d vertex_type;

	const point2d SOURCE;
	const point2d SINK;

	graphcut_graph_t(const box2d& domain, const Nbh& nbh)
	  : base(domain, nbh),
	    SOURCE ( (domain.pmin - 1) * 2 ),
	    SINK ( domain.pmax * 2 )
	{
	  m_src.resize(domain);
	  m_sink.resize(domain);
	}

	vertex_type  source(const edge_type& e) const
	{
	  switch (e.type) {
	    case SOURCE_EDGE: return SOURCE;
	    case SINK_EDGE: return e.p;
	    default: return base::source(e.p);
	  }
	}


	point2d  target(const edge_type& e) const
	{
	  switch (e.type) {
	    case SOURCE_EDGE: return e.p;
	    case SINK_EDGE: return SINK;
	    default: return base::target(e.p);
	  }
	}

	Etype& edge_at(const edge_type& e)
	{
	  if (e.type == SOURCE_EDGE)
	    return m_src.at(e.p / 2);
	  else if (e.type == SINK_EDGE)
	    return m_sink.at(e.p / 2);
	  else
	    return base::edge_at(e.p);
	}

	Etype& edge(const edge_type& e)
	{
	  mln_precondition((e.type == SOURCE_EDGE and m_src.domain().has(e.p/2)) or
			   (e.type == NORMAL_EDGE and base::m_domain.has(e.p))   or
			   (e.type == SINK_EDGE   and m_sink.domain().has(e.p/2)));
	  return edge_at(e);
	}

	Etype& edge(const point2d& v1, const point2d& v2)
	{
	  mln_precondition(v1 == SOURCE or v1 == SINK or base::_is_vertex(v1));
	  mln_precondition(v2 == SOURCE or v2 == SINK or base::_is_vertex(v2));
	  mln_precondition(v1 != v2);

	  if (v1 == SOURCE)
	    return edge(edge_type {SOURCE_EDGE, v2});
	  else if (v1 == SINK)
	    return edge(edge_type {SINK_EDGE, v2});
	  else if (v2 == SOURCE)
	    return edge(edge_type {SOURCE_EDGE, v1});
	  else if (v2 == SINK)
	    return edge(edge_type {SINK_EDGE, v1});
	  else
	    return edge(edge_type {NORMAL_EDGE, (v1 + v2) / 2});
	}


	const Etype& edge_at(const edge_type& e) const
	{
	  return const_cast<graphcut_graph_t*>(this)->edge_at(e);
	}

	const Etype& edge(const edge_type& e) const
	{
	  return const_cast<graphcut_graph_t*>(this)->edge(e);
	}

	const Etype& edge(const point2d& v1, const point2d& v2) const
	{
	  return const_cast<graphcut_graph_t*>(this)->edge(v1, v2);
	}


	Vtype& vertex(const point2d& v)
	{
	  mln_precondition(v == SOURCE or v == SINK or base::m_domain.has(v));
	  return base::vertex_at(v);
	}

	const Vtype& vertex(const point2d& v) const
	{
	  mln_precondition(v == SOURCE or v == SINK or base::m_domain.has(v));
	  return base::vertex_at(v);
	}




	vertices_range	vertices() const
	{
	  return vertices_range(vertices_iterator(base::vertices(), SOURCE, SINK));
	}

	edges_range	edges() const
	{
	  return edges_range(edges_iterator(base::edges(), base::vertices()));
	}

	adjacency_vertex_range adjacent_vertices(const point2d& v) const
	{
	  return adjacency_vertex_range(adjacency_vertex_iterator(base::adjacent_vertices(v).iter(), base::vertices(), SOURCE, SINK, v));
	}

	adjacency_edge_range	edges(const point2d& v) const
	{
	  return adjacency_edge_range(adjacency_edge_iterator(base::edges(v).iter(), base::vertices(), SOURCE, SINK, v));
	}


      private:
	image2d<Etype> m_src;
	image2d<Etype> m_sink;
      };


      // struct node_t
      // {
      // 	char lbl; // node affectation
      // 	point2d par; // parent of the node
      // 	W csource;
      // 	W csink;
      // };


      // point2d
      // findroot(const image2d<node_t>& nodes, const point2d& p)
      // {
      // 	while (nodes(p).par != p)
      // 	  p = nodes(p).par;
      // 	return p;
      // }


      template <typename G>
      bool
      is_orphan(G& graph, const point2d& p)
      {
	if (graph.vertex(p).par == p)
	  return (p != graph.SOURCE and p != graph.SINK);

	bool b = is_orphan(graph, graph.vertex(p).par);
	if (b)
	  graph.vertex(p).par = p;
	return b;
      }

      /*
      template <typename G>
      bool
      try_adopt(G& graph, const point2d& p)
      {
	auto& x = graph.vertex(p);
	mln_foreach(const point2d& q, graph.adjacet_vertices(p))
	  {
	    auto& y = graph.vertex(q);
	    bool orphan = is_orphan(graph, q);
	    if (!orphan and y.zpar == x.zpar and tree_cap(graph, p, q) > 0)
	      {
		x.par = q;
		return true;
	      }
	  }
	return false;
      }

      template <typename G>
      void depthfirst_propagate(G& graph, const point2d& x, std::queue<point2d>& actives)
      {
	std::queue<point2d> queue;
	queue.push_back(x);
	while (!queue.empty())
	  {
	    point2d p = queue.front();
	    queue.pop();
	    auto& x = graph.vertex(p);
	    mln_foreach(const point2d& q, graph.adjacent_vertices(p))
	      {
		auto& y = graph.vertex(q);
		if (y.zpar == x.zpar and tree_cap(graph, q, p) > 0 and is_orphan(graph, q))
		  {
		    y.par = p;
		    queue.push(q);
		  }
	      }
	  }
      }


      template <typename G, typename W>
      bool
      update_and_adopt(G& graph, const point2d& x, const W& delta)
      {
	point2d& y = graph.vertex(x).par;
	if (y == x) // root node
	  return true;
	else {
	  W v = (graph.edge(x, y) -= delta);
	  bool r = update_and_adopt(graph, y, delta);
	  if (v == 0)
	    y = x;

	  if (not r) // par(x) is an orphan
	    {
	      // try x to be adopted
	      if (try_adopt(graph, x))
		{
		  depthfirst_propagate(graph, x);
		}


	    }
	  if (v == 0)
	    return false;
	}
      }
      */

      template <typename DataFidelityFunction, typename RegularityFunction, typename V>
      struct default_cmp
      {
	typedef typename std::result_of<DataFidelityFunction(bool, V)>::type W1;
	typedef typename std::result_of<RegularityFunction(bool, bool)>::type W2;

	typedef typename std::common_type<W1, W2>::type W;

	int
	operator () (W x, W y) const
	{
	  return (y < x) - (x < y);
	}
      };

    }


    template <typename V, class N, class DataFidelityFunction, class RegularityFunction, class Compare>
    double
    graphcut(const image2d<V>& ima, image2d<bool>& f, const N& nbh, DataFidelityFunction d, RegularityFunction v, const Compare& cmp)
    {
      typedef typename std::result_of<DataFidelityFunction(bool, V)>::type W1;
      typedef typename std::result_of<RegularityFunction(bool, bool)>::type W2;

      typedef typename std::common_type<W1, W2>::type W;


      struct node_t
      {
	point2d par;
	point2d zpar;
	bool	active; // active means also INQUEUE
      };

      // a node x is free iif x != SINK, SOURCE, x.par = x, x.zpar = x , not(x.active)
      // a node x is active iif x.zpar = (SOURCE|SINK) and x.active
      // a node x is an orphan iif x.par = x and x.zpar = (SOURCE|SINK)

      typedef internal::graphcut_graph_t<node_t, std::pair<W,W>, N>  G;

      trace::entering("mln::graphcut::graphcut");

      G graph(ima.domain(), nbh);

      // initialize graph
      {
	// Constraints on p -> q
	W w,w2 = v(true, false);
	mln_foreach(auto e, graph.edges())
	  switch (e.type) {
	    case G::SOURCE_EDGE:
	      w = d(false, ima(graph.target(e) / 2));
	      graph.edge(e) = { w, 0 };
	      break;
	    case G::NORMAL_EDGE:
	      graph.edge(e) = { w2, w2 };
	      break;
	    case G::SINK_EDGE:
	      w = d(true, ima(graph.source(e) / 2));
	      graph.edge(e) = { w, 0 };
	      break;
	  }

	mln_foreach(const point2d& p, graph.vertices())
	  graph.vertex(p) = {p, p, false};
      }

      std::queue<point2d> active;
      active.push(graph.SOURCE);
      active.push(graph.SINK);
      graph.vertex(graph.SOURCE).active = true;
      graph.vertex(graph.SINK).active = true;

      auto getedge = [&graph] (const point2d& p, const point2d& q) -> W& {
	return (p < q) ? graph.edge(p,q).first : graph.edge(p,q).second;
      };

      auto tree_cap = [&graph, &getedge] (const point2d& p, const point2d& q) -> W& {
	point2d r = graph.vertex(p).zpar;
	mln_assertion(r == graph.SOURCE or r == graph.SINK);
	return (r == graph.SOURCE) ? getedge(p, q) : getedge(q, p);
      };

      double maxflow = 0;
      std::vector<point2d> orphans;

      point2d oldp, p, zp, q, zq;
      oldp = graph.SOURCE;
      p = graph.SOURCE;
      zp = graph.SOURCE;
      mln_iter(qit, graph.adjacent_vertices(p));
      qit.init();
      while (true)
	{
	  // grow step: find an augmenting path P from s to t
	  while (! active.empty())
	    {
	      p = oldp;
	      if (active.front() != p)
		{
		  //std::cout << "oldp: " << p << std::endl;
		  p = active.front();
		  zp = graph.vertex(p).zpar;
		  //std::cout << "Get " << p << "from active list. (zpar=" << zp << ")" << std::endl;
		  qit.init();
		}
	      else {
		//std::cout << "Continue " << p << "from active list. (zpar=" << zp << ")" << std::endl;
	      }
	      oldp = p;
	      if (not graph.vertex(p).active) {
		//std::cout << "  Dismissed." << std::endl;
		active.pop();
		continue;
	      }

	      for(;not qit.finished(); qit.next()) {
		q = *qit;
		if (cmp(tree_cap(p, q), 0) > 0)
		    {
		      zq = graph.vertex(q).zpar;
		      if (zq == q and zq != graph.SOURCE and zq != graph.SINK) // q is a free node, add it as an an active node
			{
			  bool inqueue = graph.vertex(q).active;
			  graph.vertex(q) = {p, zp, true};
			  if (!inqueue)
			    active.push(q);
			  //std::cout << "  Push " << q << "(par,zpar,cap: " << p << "," << zp << "," << tree_cap(p,q) << ") in active list." << std::endl;
			}
		      else if (zq != zp) // we have found a path from SOURCE to SINK
			goto augmentation;
		    }
	      }
	      // p is not an active node anymore
	      graph.vertex(p).active = false;
	      active.pop();
	    }

	  // No path has been found
	  break;

	augmentation:
	  //std::cout << "Found path: " << p << " -> " << q << std::endl;
	  mln_assertion(graph.vertex(p).zpar == graph.SOURCE or graph.vertex(p).zpar == graph.SINK);
	  mln_assertion(graph.vertex(q).zpar == graph.SOURCE or graph.vertex(q).zpar == graph.SINK);
	  mln_assertion(graph.vertex(p).zpar != graph.vertex(q).zpar);

	  if (graph.vertex(q).zpar == graph.SOURCE)
	    std::swap(p,q);

	  // find bottleneck
	  W delta = value_traits<W>::max();
	  {
	    for (point2d x: {p, q})
	      {
		point2d y = graph.vertex(x).par;
		while (x != y) {
		  //std::cout << x << "-->" << y << std::endl;
		  delta = std::min(delta, tree_cap(y, x));
		  x = y;
		  y = graph.vertex(x).par;
		}
	      }
	    delta = std::min(delta, tree_cap(p, q));
	  }

	  //std::cout << "  -- Delta: " << delta << std::endl;
	  maxflow += delta;

	  // Update path
	  {
	    orphans.clear();
	    {
	      point2d x = p;
	      point2d y = graph.vertex(x).par;
	      while (x != y) {
		//std::cout << "  up(" << x << "," << y << "):" << getedge(y,x) << std::endl;
		W v = (getedge(y, x) -= delta);
		getedge(x, y) += delta;
		  //std::cout << "  up(" << x << "," << y << "):" << getedge(y,x) << std::endl;
		if (cmp(v, 0) == 0) { // x->y become saturated
		  orphans.push_back(x);
		  graph.vertex(x).par = x;
		}
		x = y;
		y = graph.vertex(x).par;
	      }
	    }
	    {
	      point2d x = q;
	      point2d y = graph.vertex(x).par;
	      while (x != y) {
		  W v = (getedge(x, y) -= delta);
		  getedge(y, x) += delta;
		  if (cmp(v, 0) == 0) { // x->y become saturated
		    orphans.push_back(x);
		    graph.vertex(x).par = x;
		  }
		  x = y;
		  y = graph.vertex(x).par;
	      }
	    }

	    W v = (getedge(p, q) -= delta);
	    getedge(q, p) += delta;
	    //std::cout << "e(" << p << "," << q << "):" << getedge(p,q) << std::endl;
	  }


	  // Adoption stage
	  {
	    point2d p = {0,0}; 
	    mln_iter(qit, graph.adjacent_vertices(p));

	    for (unsigned i = 0; i < orphans.size(); ++i)
	      {
		p = orphans[i];
		node_t& x = graph.vertex(p);
		//std::cout << " -- Proc orphan " << p << "(par=" << x.par << ",zpar=" << x.zpar << ")" << std::endl;
		bool adopted = false;
		mln_forall(qit)
		  {
		    const point2d& q = *qit;
		    bool orphan = internal::is_orphan(graph, q);
		    node_t& y = graph.vertex(q);
		    if (!orphan and x.zpar == y.zpar and cmp(tree_cap(q, p), 0) > 0)
		      {
			//std::cout << "    adopted by" << q << "(par=" << y.par << ",zpar=" << y.zpar << ")" << std::endl;
			x.par = q; // y adopt x
			adopted = true;
			break;
		      }
		  }
		if (!adopted)
		  {
		    mln_forall(qit)
		      {
			const point2d& q = *qit;
			node_t& y = graph.vertex(q);
			if (x.zpar == y.zpar) {
			  if (cmp(tree_cap(q, p), 0) > 0) { // p's neighbors become active
			    if (!y.active) {
			      active.push(q);
			      y.active = true;
			    }
			  }
			  if (y.par == p) { // q points to an orphan so is an orphan itself
			    //std::cout << "    Setting oprhan " << q << std::endl;
			    y.par = q;
			    orphans.push_back(q);
			  }
			}
		      }

		    // p becomes a free node
		    x.par = p;
		    x.zpar = p;
		    x.active = false;
		  }
	      }
	  }
	  //std::cout << orphans.size() << std::endl;
	}


      mln_pixter(px, f);
      mln_forall(px)
      {
	point2d p = graph.vertex(px->point() * 2).zpar;
	//mln_assertion(p == graph.SOURCE or p == graph.SINK);
	px->val() = (p != graph.SINK);
      }

      std::cout << "Maxflow: " << maxflow << std::endl;

      trace::exiting();
      return maxflow;
    }

  }

}

#endif // ! GRAPHCUT_HH
