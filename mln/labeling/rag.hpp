#ifndef MLN_LABELING_RAG_HPP
# define MLN_LABELING_RAG_HPP

# include <boost/graph/adjacency_list.hpp>
# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/extension/extension.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/core/trace.hpp>
# include <utility>
# include <vector>

namespace mln
{
  namespace labeling
  {

    struct vertex_label_t
    {
      typedef boost::vertex_property_tag kind;
    };

    namespace internal
    {
      template <typename Label>
      struct rag_graph
      {
	typedef boost::property<vertex_label_t, Label> Vproperty;
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vproperty> type;
      };

    }

    template <typename I, typename N, typename Graph, typename Label = unsigned>
    std::pair<mln_ch_value(I, Label), Label>
    rag(const Image<I>& ima, const N& nbh, Graph& graph, Label lbl = unsigned());


/******************************************/
/****          Implementation          ****/
/******************************************/

    namespace impl
    {


      template <typename I, typename N, typename Graph, typename Label, typename J>
      Label
      rag(const I& ima, const N& nbh, Graph& graph, Label lbl, J&& out)
      {
        typedef mln_value(I) V;
        typedef mln_point(I) P;

        const Label bg = value_traits<Label>::max();
        const Label INFRONT = lbl;

        extension::fill(out, INFRONT);

        P* queue = new P[ima.domain().size()];

        unsigned qstart = 0, qend = 0;
        unsigned fstart = ima.domain().size();
        unsigned fend = ima.domain().size();

        auto qpush = [&queue, &qstart, &qend] (const P& p) { queue[qend++] = p; };
        auto qpop = [&queue, &qstart, &qend] () -> P { return queue[--qend]; };
        auto qempty = [&queue, &qstart, &qend] () -> bool { return qend == 0; };

        auto fpush = [&queue, &fstart, &fend] (const P& p) { queue[--fstart] = p; };
        auto fpop =  [&queue, &fstart, &fend] () -> P { return queue[--fend]; };
        auto fempty = [&queue, &fstart, &fend] () -> bool { return fend == fstart; };

        P q;
        mln_iter(n, nbh(q));

        mln_foreach(P p, ima.domain())
        {
          if (ima(p) and out(p) == bg)
            {
              qpush(p);
              out(p) = ++lbl;
              mln_assertion(lbl < value_traits<Label>::max());
              while (not qempty())
                {
                  q = qpop();
                  mln_forall(n)
                    if (out.at(*n) == bg)
                      {
                        if (ima(*n)) {
                          qpush(*n);
                          out(*n) = lbl;
                        } else {
                          fpush(*n);
                          out(*n) = INFRONT;
                        }
                      }
                }
            }
        }

	graph = Graph(lbl);

	// Add vertices to the graph
	typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
	std::vector< typename boost::graph_traits<Graph>::vertex_descriptor > lbl2vertex;
	lbl2vertex.resize(lbl+1 - INFRONT);

	std::tie(vi, vi_end) = boost::vertices(graph);
	for (Label i = INFRONT + 1; vi != vi_end; ++vi, ++i) {
	  boost::put(vertex_label_t(), graph, *vi, i);
	  lbl2vertex[i - INFRONT] = *vi;
	}

	std::vector<Label> tomerge;
        while (not fempty())
          {
            q = fpop();
            Label& v = out(q);
	    bool removed = false;
	    tomerge.clear();
            mln_forall(n)
            {
              Label x = out.at(*n);
              if (x != INFRONT and x != bg)
                {
                  if (!removed and v == INFRONT)
                    v = x;
                  else if (removed or v != x) {
		    removed = true;
		    tomerge.push_back(x);
                  }
                }
            }

            if (removed) {
	      tomerge.push_back(v);
	      int n = tomerge.size();
	      for (int i = 0; i < n; ++i)
		for (int j = i+1; j < n; ++j)
		  if (tomerge[i] != tomerge[j])
		    boost::add_edge(lbl2vertex[tomerge[i] - INFRONT], lbl2vertex[tomerge[j] - INFRONT], graph);
	      v = INFRONT;
              continue;
            }

            mln_forall(n)
              if (out.at(*n) == bg)
                {
                  fpush(*n);
                  out(*n) = INFRONT;
                }
          }

        mln_pixter(px, out);
        mln_forall(px)
          if (px->val() == bg)
            px->val() = INFRONT;


        delete [] queue;
        return lbl;
      }


    } // end of namespace mln::labeling::impl


    template <typename I, typename N, typename Graph, typename Label>
    std::pair<mln_ch_value(I, Label), Label>
    rag(const Image<I>& ima_, const N& nbh, Graph& graph, Label lbl)
    {
      static_assert(std::is_convertible<mln_value(I), bool>::value,
                    "Image value type must be convertible to bool.");
      static_assert(std::is_same<typename boost::graph_traits<Graph>::directed_category, boost::undirected_tag>::value,
		    "Graph must undirected");
      static_assert(std::is_same<typename boost::graph_traits<Graph>::edge_parallel_category, boost::disallow_parallel_edge_tag>::value,
		    "Graph must not allow parallel edges.");
      static_assert(boost::has_vertex_property<Graph>::value, "Graph must have a vertex property map");

      typedef typename boost::property_map<Graph, vertex_label_t>::type PM;
      static_assert(std::is_convertible<Label, typename boost::property_traits<PM>::value_type>::value,
		    "Graph must have an internal vertex label property map and should be convertible to Label");

      mln_entering("mln::labeling::rag");

      const I& ima = exact(ima_);
      const Label bg = value_traits<Label>::max();

      int status;
      mln_ch_value(I, Label) out = imchvalue<Label>(ima)
        .adjust(nbh)
        .init(bg)
        .get_status(status);

      if (status == 0)
        lbl = impl::rag(ima, nbh, graph, lbl, out);
      else
        lbl = impl::rag(ima, nbh, graph, lbl, extension::add_value_extension(out,bg));

      mln_exiting();

      return std::make_pair(out, lbl);
    }

  } // end of namespace mln::labeling

} // end of namespace mln

#endif //!MLN_LABELING_RAG_HPP
