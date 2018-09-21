#ifndef SIMPLIFY_HPP
# define SIMPLIFY_HPP

#include <vector>

#include <mln/core/image/image2d.hpp>
#include <mln/core/extension/fill.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/win2d.hpp>
#include <apps/tos/topology.hpp>


namespace mln
{




  /// \brief Simplify an image removing overlapping level lines.
  ///
  /// \param f: The image on which the tree has been computed
  /// \param K: K outputed by the ToS algorithm (twice as big as \p f)
  /// \param parent: parent outputed by the ToS algorithm
  /// \param S: S outputed by the ToS algorithm
  /// \param lambda: radius for non overlapping level lines
  /// \param area: perform a area grain filter
  /// \param areafactor: avoid sur-simplification by limiting node desactivation
  ///                    a shape S1 can only be desactived by a shape S2 if area(S1) * areafactor < area(S2)
  template <typename V, typename T>
  image2d<V>
  simplify_bottom_up(const image2d<V>& f, const image2d<T>& K,
		     const image2d<unsigned>& parent, const std::vector<unsigned>& S,
		     int lambda, int area = 0, float areafactor = 0);


  template <typename V, typename T>
  image2d<V>
  simplify_top_down(const image2d<V>& f, const image2d<T>& K,
		    const image2d<unsigned>& parent, const std::vector<unsigned>& S, int eps);

  /// \brief For any shape \p S, retrieve the smallest shape \p A such that $\delta_\epslison(s) \subset A$
  ///
  /// A stucturing element used for dilation is square of radius $\epsilon$.
  template <typename T>
  image2d<unsigned>
  smallest_enclosing_shape(const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S, int eps,
			   const image2d<unsigned>* depth = NULL);


  template <typename T>
  image2d<unsigned>
  compute_depth(const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S);


  /*******************************/
  /** Implementation            **/
  /*******************************/

  namespace internal
  {

    template <typename T, typename N>
    unsigned
    common_ancestor(unsigned x,
		    const Neighborhood<N>& nbh_,
		    const image2d<T>& K,
		    const image2d<unsigned>& parent,
		    const image2d<unsigned>& depth)
    {
      std::vector<unsigned> v;
      unsigned minp = K[x] == K[parent[x]] ? parent[x] : x;
      unsigned mindepth = depth[minp];
      v.push_back(minp);


      const N& nbh = exact(nbh_);
      point2d p = K.point_at_index(x);
      mln_foreach(const point2d& n, nbh(p))
	{
	  if (K.domain().has(n))
	    {
	      unsigned q = K.index_of_point(n);

	      // WARNING: Le root est pute, dans le cas du thickening
	      // il ne faut pas le considérer car de nombreuses ont au moins un edge
	      // en commun avec lui (border effect).
	      unsigned rq = K[q] == K[parent[q]] ? parent[q] : q;
 	      // if (parent[rq] == rq) // Root -> casse toi
	      // 	continue;

	      v.push_back(rq);
	      if (depth[rq] < mindepth) {
		mindepth = depth[rq];
		minp = rq;
	      }
	    }
	}

      bool modif;
      do {
	modif = false;
	for(unsigned& x: v)
	  if (depth[x] > mindepth) {
	    x = parent[x];
	    modif = true;
	  } else if (depth[x] == mindepth and x != minp) {
	    x = parent[x];
	    mindepth--;
	    minp = x;
	    modif = true;
	  }
      } while (modif);

      mln_assertion(std::all_of(v.begin(), v.end(), [minp](unsigned x) { return x == minp; }));
      return minp;
    }



    unsigned
    common_ancestor(unsigned x, unsigned y, const image2d<unsigned>& parent, const image2d<unsigned>& depth)
    {
      while (x != y) {
	if (depth[x] > depth[y])
	  x = parent[x];
	else
	  y = parent[y];
      }
      return x;
    }
  }


  template <typename T>
  image2d<unsigned>
  compute_depth(const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    image2d<unsigned> depth;
    resize(depth, K);

    // Compute depth attribute
    {
      depth[S[0]] = 0;
      for (unsigned i = 1; i < S.size(); ++i)
	{
	  unsigned x = S[i];
	  if (K[x] != K[parent[x]]) // canonical element
	    depth[x] = depth[parent[x]] + 1;
	  else
	    depth[x] = depth[parent[x]];
	}
    }
    return depth;
  }



  template <typename T>
  image2d<unsigned>
  smallest_enclosing_shape(const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S, int eps,
			   const image2d<unsigned>* depth_)
  {

    image2d<unsigned> depth = (depth_ != NULL) ? (*depth_) : compute_depth(K, parent, S);

    image2d<unsigned> enc;
    resize(enc, parent);


    // Init at beginning, the most enclosing shape of S is S itself.
    mln_pixter(px, enc);
    mln_forall(px) {
      unsigned i = px->index();
      px->val() = (K[i] == K[parent[i]]) ? parent[i] : i;
    }

    // Nbh
    rect2d win = make_rectangle2d(2*eps+1, 2*eps+1);

    for (int i = S.size()-1; i > 0; --i)
      {
	unsigned x = S[i];

	if (not K1::is_face_2(K.point_at_index(x)))
	  continue;

	unsigned anc = internal::common_ancestor(x, win, K, parent, depth);
	enc[x] = internal::common_ancestor(enc[x], anc, parent, depth);
	enc[parent[x]] = internal::common_ancestor(enc[parent[x]], enc[x], parent, depth);
      }

    return enc;
  }




  template <typename V, typename T>
  image2d<V>
  simplify_bottom_up(const image2d<V>& f,
		     const image2d<T>& K,
		     const image2d<unsigned>& parent,
		     const std::vector<unsigned>& S,
		     int eps,
		     int area,
		     float areafactor)
  {
    image2d<unsigned> areas;
    resize(areas, K).init(0);

    // Compute node activity
    image2d<bool> active;
    resize(active, K).init(true);


    // Compute area
    for (int i = S.size() - 1; i > 0; --i)
      {
	unsigned x = S[i];
	if (not K1::is_face_2(K.point_at_index(x)))
	  continue;
	areas[x] += 1;
	areas[parent[x]] += areas[x];
	active[x] = areas[x] >= area;
      }
    // root
    areas[S[0]] += 1;

    image2d<unsigned> enc = smallest_enclosing_shape(K, parent, S, eps);
    image2d<unsigned> depth;
    resize(depth, K);
    // Compute depth attribute
    {
      depth[S[0]] = 0;
      for (unsigned i = 1; i < S.size(); ++i)
	{
	  unsigned x = S[i];
	  if (K[x] != K[parent[x]]) // canonical element
	    depth[x] = depth[parent[x]] + 1;
	  else
	    depth[x] = depth[parent[x]];
	}
    }


    int nactive = 0;
    for (int i = S.size() - 1; i > 0; --i)
      {
	unsigned x = S[i];

	if (K[x] == K[parent[x]] or !active[x])
	  continue;

	++nactive;

	unsigned anc = enc[x];
	//std::cout << "Node depth: " << depth[x] << " -> Ancestor depth: " << depth[anc] << std::endl;

	// unactive ]x ---> anc[
	if (x != anc) {
	  unsigned p = parent[x];
	  while (p != anc and areas[p] * areafactor < areas[x])
	    {
	      active[p] = false;
	      p = parent[p];
	    }
	}
      }

    // Simplify
    image2d<V> out;
    resize(out, f);
    nactive = 0;
    for (unsigned x: S)
      {
	point2d p = K.point_at_index(x);
	if (K1::is_face_2(p))
	  {
	    if ((active[x] and K[parent[x]] != K[x]) or x == parent[x]) {
	      ++nactive;
	      out(p/2) = f(p/2);
	    } else {
	      point2d q = K.point_at_index(parent[x]);
	      mln_assertion(K1::is_face_2(q));
	      out(p/2) = out(q/2);
	    }
	  }
      }
    std::cout << "Number of nodes: " << nactive << std::endl;
    //io::imsave(transform(enc, [&areas, &S](unsigned x) -> float { return (float) areas[x] / areas[S[0]]; }), "tmp.tiff");

    return out;
  }

  template <typename V, typename T, class Predicate>
  image2d<V>
  simplify_top_down_pred(const image2d<V>& f, const image2d<T>& K, const image2d<unsigned>& parent,
			 const std::vector<unsigned>& S, Predicate pred)
  {
    // Compute node activity
    image2d<unsigned> alive;
    resize(alive, K);

    alive[S[0]] = S[0];

    for (int i = 1; i < S.size(); ++i)
      {
	unsigned x = S[i];
	if (K[x] == K[parent[x]])
	  continue;

	unsigned y = parent[x];
	unsigned z = alive[y]; // last ancestor alive
	alive[x] = pred(x, z) ? x : z;
      }

    // Simplify
    image2d<V> out;
    resize(out, f);
    int nactive = 0;

    // point2d p = K.point_at_index(S[0]);
    // out(p/2) = f(p/2);
    // for (unsigned x: S)
    //   {
    // 	point2d p = K.point_at_index(x);
    // 	if (K1::is_face_2(p))
    // 	  {
    // 	    if ((active[x] and K[parent[x]] != K[x]) or x == parent[x]) {
    // 	      ++nactive;
    // 	      out(p/2) = f(p/2);
    // 	    } else {
    // 	      point2d q = K.point_at_index(parent[x]);
    // 	      mln_assertion(K1::is_face_2(q));
    // 	      out(p/2) = out(q/2);
    // 	    }
    // 	  }
    //   }
    // std::cout << "Number of nodes: " << nactive << std::endl;
    // std::cout << "Number of nodes with root connection:" << root_connection << std::endl;
    //io::imsave(transform(enc, [&areas, &S](unsigned x) -> float { return (float) areas[x] / areas[S[0]]; }), "tmp.tiff");

    return out;


  }


  template <typename V, typename T>
  image2d<V>
  simplify_top_down_tolerance(const image2d<V>& f, const image2d<T>& K,
			      image2d<unsigned>& parent,
			      const std::vector<unsigned>& S, int eps, float tol)
  {
    // Smallest enclosing shape
    image2d<unsigned> depth = compute_depth(K, parent, S);
    image2d<unsigned> enc = smallest_enclosing_shape(K, parent, S, eps, &depth);

    // Attribute
    struct attr_t
    {
      std::set<point2d> ext;
      std::vector<unsigned> inter;


      std::ostream&
      pprint (std::ostream& os) {
	os << "Inter / Ext = [ ";
	for (unsigned v: inter)
	  os << v << " ";
	return os << " ] / " << ext.size();
      }
    };

    // Set a fake root
    unsigned ROOT = (unsigned) -1;
    parent[S[0]] = ROOT;

    image2d<attr_t> attr;
    resize(attr, K);
    {
      auto realnode = [&K, &parent, ROOT](unsigned x) {
	return (parent[x] == ROOT or K[x] != K[parent[x]]) ? x : parent[x];
      };

      auto se = make_rectangle2d(2*eps+1, 2*eps+1);

      mln_iter(p_, K.domain());
      mln_iter(n_, se(p_));
      mln_foreach(const point2d& p, p_)
      {
	unsigned pnode = parent.index_of_point(p);
	pnode = realnode(pnode);

	mln_foreach(const point2d& n, n_)
	{
	  unsigned qnode, nnode;
	  //std::cout << n << std::endl;

	  if (parent.domain().has(n))
	    {
	      nnode = parent.index_of_point(n);
	      nnode = realnode(nnode);
	      qnode = internal::common_ancestor(pnode, nnode, parent, depth);
	    }
	  else
	    {
	      nnode = ROOT;
	      qnode = ROOT;
	    }

	  // n ∈ Ext(S, ε) with S ∈ [Sp  → Sq[
	  // n ∈ S' with S' ∈ [Sq  → root [ (but ENCLOSING(S) is enough)
	  // ∀ (S,S') n ∈ Ext(S, ε) ⋂ S' ≠ ∅
	  unsigned x = pnode;
	  bool inserted = false;
	  while (x != qnode and !inserted)
	    {
	      std::tie(std::ignore, inserted) = attr[x].ext.insert(n);
	      if (inserted)
		{
		  //std::cout << "insert: " << n << " in ext-" << x << std::endl;
		  unsigned y = parent[x];
		  unsigned i = 0;
		  while (y != qnode) {
		    ++i;
		    y = parent[y];
		  }
		  if (attr[x].inter.size() < i)
		    attr[x].inter.resize(i);

		  while (y != enc[x] and y != ROOT)
		    {
		      if (i < attr[x].inter.size())
			++attr[x].inter[i];
		      else
			attr[x].inter.push_back(1);
		      ++i;
		      y = parent[y];
		    }
		}
	      x = parent[x];
	    }
	}
      }
    }

    // Compute node activity
    image2d<bool> active;
    resize(active, K).init(false);

    active[S[0]] = true;
    //std::cout << "Activate: " << S[0] << std::endl;
    auto crit = [&attr, tol](unsigned p, unsigned k) {
      return k >= attr[p].inter.size() or (attr[p].inter[k] > (tol * attr[p].ext.size()));
    };

    for (unsigned i = 1; i < S.size(); ++i)
      {
	unsigned x = S[i];
	if (K[x] == K[parent[x]])
	  continue;

	unsigned anc = enc[x];
	// Prevent bordering effect (if anc == root, set x active)
	// if (anc == S[0]) {
	//   root_connection++;
	//   active[x] = true;
	//   continue;
	// }

	unsigned p = parent[x];
	unsigned k = 0;
	// A node is active if ]x --> y[ are inactive and y
	// verfies the criterion
	while (p != anc and !active[p] and !crit(x, k)) {
	  p = parent[p];
	  ++k;
	}

	//std::cout << x << " --> " << parent[x] << " / k=" << k << " !";
	//attr[k].pprint(std::cout) << " res=" << crit(x,k) << std::endl;

	if (crit(x,k)) {
	  // std::cout << "Activate: " << x << " in " << p
	  // 	    << " Anc=" << anc << " K=" << k << std::endl;
	  active[x] = true;
	}
      }

    // Restore root
    parent[S[0]] = S[0];

    // Simplify
    image2d<V> out;
    resize(out, f);
    int nactive = 0;
    int nbefore = 0;

    point2d p = K.point_at_index(S[0]);
    out(p/2) = f(p/2);
    for (unsigned x: S)
      {
	if (K[parent[x]] != K[x] or x == parent[x])
	  nbefore++;

	point2d p = K.point_at_index(x);
	if (K1::is_face_2(p))
	  {
	    if ((active[x] and K[parent[x]] != K[x]) or x == parent[x]) {
	      ++nactive;
	      out(p/2) = f(p/2);
	    } else {
	      point2d q = K.point_at_index(parent[x]);
	      mln_assertion(K1::is_face_2(q));
	      out(p/2) = out(q/2);
	    }
	  }
      }
    std::cout << "Number of nodes (before): " << nbefore << std::endl;
    std::cout << "Number of nodes (after): " << nactive << std::endl;
    //std::cout << "Number of nodes with root connection:" << root_connection << std::endl;
    //io::imsave(transform(enc, [&areas, &S](unsigned x) -> float { return (float) areas[x] / areas[S[0]]; }), "tmp.tiff");

    return out;

  }



  template <typename V, typename T>
  image2d<V>
  simplify_top_down(const image2d<V>& f, const image2d<T>& K,
   		    const image2d<unsigned>& parent, const std::vector<unsigned>& S, int eps)
  {
    // Compute node activity
    image2d<bool> active;
    resize(active, K).init(false);


    image2d<unsigned> enc = smallest_enclosing_shape(K, parent, S, eps);
    active[S[0]] = true;
    //std::cout << "Activate: " << S[0] << std::endl;
    unsigned root_connection = 0;
    for (unsigned i = 1; i < S.size(); ++i)
      {
	unsigned x = S[i];
	if (K[x] == K[parent[x]])
	  continue;

	unsigned anc = enc[x];
	// Prevent bordering effect (if anc == root, set x active)
	if (anc == S[0]) {
	  root_connection++;
	  // active[x] = true;
	  // continue;
	}

	unsigned p = x;
	if (p != anc) {
	  // A node is active if ]x --> anc[ are inactive
	  p = parent[x];
	  while (p != anc and !active[p])
	    p = parent[p];
	}

	if (p == anc) {
	  //std::cout << "Activate: " << x << " in " << p << " Anc= " << anc << std::endl;
	  active[x] = true;
	}
      }

    // Simplify
    image2d<V> out;
    resize(out, f);
    int nactive = 0;

    point2d p = K.point_at_index(S[0]);
    out(p/2) = f(p/2);
    for (unsigned x: S)
      {
	point2d p = K.point_at_index(x);
	if (K1::is_face_2(p))
	  {
	    if ((active[x] and K[parent[x]] != K[x]) or x == parent[x]) {
	      ++nactive;
	      out(p/2) = f(p/2);
	    } else {
	      point2d q = K.point_at_index(parent[x]);
	      mln_assertion(K1::is_face_2(q));
	      out(p/2) = out(q/2);
	    }
	  }
      }
    std::cout << "Number of nodes: " << nactive << std::endl;
    std::cout << "Number of nodes with root connection:" << root_connection << std::endl;
    //io::imsave(transform(enc, [&areas, &S](unsigned x) -> float { return (float) areas[x] / areas[S[0]]; }), "tmp.tiff");

    return out;

  }

}

#endif // ! SIMPLIFY_HPP
