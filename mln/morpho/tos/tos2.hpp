// Copyright (C) 2009 EPITA Research and Development Laboratory (LRDE)
//
// This file is part of Olena.
//
// Olena is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, version 2 of the License.
//
// Olena is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Olena.  If not, see <http://www.gnu.org/licenses/>.
//
// As a special exception, you may use this file as part of a free
// software project without restriction.  Specifically, if other files
// instantiate templates or use macros or inline functions from this
// file, or you compile this file and link it with other files to produce
// an executable, this file does not by itself cause the resulting
// executable to be covered by the GNU General Public License.  This
// exception does not however invalidate any other reasons why the
// executable file might be covered by the GNU General Public License.

#ifndef MLN_MORPHO_TOS_TOS_HPP
# define MLN_MORPHO_TOS_TOS_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/value/value_traits.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/extension/fill.hpp>
# include <mln/morpho/tos/irange.hpp>
# include <mln/morpho/tos/immerse.hpp>
# include <mln/morpho/tos/pset.hpp>


namespace mln
{

  namespace morpho
  {

    namespace internal
    {
      template <typename Compare>
      struct equiv;
    }


    /// \brief Compute the tree of shapes
    ///
    /// Compute the tree of shapes my immerging the image in Khalimsky grid.
    /// The user is responsible to provide an image with a constant border, since the
    /// algorithm starts flooding at point ima.domain().pmin. The user may also want
    /// to perform a median interpolation to ensure correct properties on TOS.
    ///
    /// \return
    ///
    template <typename I, typename Domain, typename Neighborhood, typename Compare = std::less<mln_value(I)>, typename Equiv = internal::equiv<Compare> >
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    subToS(const Image<I>& ima, const Domain& subdomain, const Neighborhood& nbh, const Compare& cmp, const Equiv& eq);

    template <typename I, typename Domain, typename Neighborhood, typename Compare = std::less<mln_value(I)> >
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    subToS(const Image<I>& ima, const Domain& subdomain, const Neighborhood& nbh, const Compare& cmp = Compare () );




    /********************/
    /** Implementation **/
    /********************/

    namespace internal
    {
      template <typename I>
      inline
      typename I::size_type
      zfind_root(I& zpar, typename I::size_type x)
      {
	if (zpar[x] != x)
	  zpar[x] = zfind_root(zpar, zpar[x]);
	return zpar[x];
      }


      template <typename Compare>
      struct equiv
      {
	equiv(const Compare& cmp) :
	m_cmp (cmp)
	{
	}

	template <typename T>
	bool operator () (const T& x, const T& y) const
	{
	  return !m_cmp(x,y) and !m_cmp(y,x);
	}

      private:
	Compare m_cmp;
      };
    }


    template <typename I, typename Domain, typename Neighborhood, typename Compare = std::less<mln_value(I)>, typename Equiv = internal::equiv<Compare> >
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    subToS(const Image<I>& ima_,const image2d<bool>& mask, const Domain& subdomain, const Neighborhood& nbh, const Compare& cmp, const Equiv& eq)
    {
      using namespace mln::morpho::tos;

      typedef mln_value(I) V;
      typedef irange<V> R;
      typedef typename I::size_type size_type;

      static constexpr size_type UNPROCESSED = value_traits<size_type>::max();
      static constexpr size_type PROCESSED = 0;

      const I& ima = exact(ima_);


      // f: image of interval in Khalimsky space
      // K: image of value in Khalimsky that tells at which a level a point is inserted
      mln_ch_value(I, R) f = tos::internal::immerse(ima, cmp);
      mln_concrete(I) K;
      mln_ch_value(I, size_type) parent, zpar;
      std::vector<size_type> S;

      //S.reserve(f.domain().size());
      resize(K, f);
      resize(parent, f).init(PROCESSED);
      extension::fill(parent, PROCESSED);

      mln_foreach(point2d p, subdomain)
	{
	  point2d q = 2*p;
	  mln_iter(n1, c8(p));
	  mln_iter(n2, c8(q)); // Q: why 2*p does not work ? A: because not a lvalue

	  mln_forall(n1, n2)
	    if (mask.at(*n1))
	      parent.at(*n2) = UNPROCESSED;
	  parent(q) = UNPROCESSED;
	}

      pset<I, Compare> W(K, cmp);
      auto pit = subdomain.iter(); pit.init();
      auto pmin = 2 * (*pit);
      size_type p = f.index_of_point(pmin);
      W.insert(p);
      parent[p] = PROCESSED;
      K[p] = f[p].lower;

      auto dindexes = wrt_delta_index(f, nbh.dpoints);
      while (!W.empty())
	{
	  //std::cout << "Search: " << K[p] << std::endl;
	  p = W.has_next(p) ? W.pop_next(p) : W.pop_previous(p);
	  //std::cout << "Found: " << p << " @ " << K[p] << std::endl;
	  V curlevel = K[p];
	  S.push_back(p);

	  mln_foreach (int k, dindexes)
	    {
	      size_type q = p + k;
	      if (parent[q] == UNPROCESSED)
		{
		  if (cmp(f[q].upper, curlevel))
		    K[q] = f[q].upper;
		  else if (cmp(curlevel, f[q].lower))
		    K[q] = f[q].lower;
		  else
		    K[q] = curlevel;

		  parent[q] = PROCESSED;
		  //std::cout << "Insert:" << q << " @ " << K[q] << std::endl;
		  W.insert(q);
		}
	    }
	}

      // 2nd step: union-find

      resize(zpar, parent).init(UNPROCESSED);
      extension::fill(zpar, UNPROCESSED);
      for (int i = S.size()-1; i >= 0; --i)
	{
	  size_type p = S[i];
	  parent[p] = p;
	  zpar[p] = p;

	  mln_foreach (int k, dindexes)
	    {
	      size_type q = p + k;
	      if (zpar[q] != UNPROCESSED)
		{
		  size_type r = morpho::internal::zfind_root(zpar, q);
		  if (r != p) {
		    parent[r] = p;
		    zpar[r] = p;
		  }
		}
	    }

	}

      // 3rd step: canonicalization
      for (size_type p : S)
	{
	  size_type q = parent[p];
	  if (eq(K[q], K[parent[q]]))
	    parent[p] = parent[q];
	}


      //mln_postcondition(S.size() == K.domain().size());
      //mln_postcondition(S.size() == parent.domain().size());
      // All done !
      //std::cout << "Size: " << S.size() << std::endl;
      return std::make_tuple(std::move(K), std::move(parent), std::move(S));
    }


    template <typename I, typename Domain, typename Neighborhood, typename Compare = std::less<mln_value(I)> >
    std::tuple< mln_concrete(I), mln_ch_value(I, typename I::size_type), std::vector<typename I::size_type> >
    subToS(const Image<I>& ima, const image2d<bool>& mask, const Domain& subdomain, const Neighborhood& nbh, const Compare& cmp = Compare () )
    {
      return subToS(ima, mask, subdomain, nbh, cmp, internal::equiv<Compare> (cmp));
    }

  }

}

#endif // ! MLN_MORPHO_TOS_TOS_HPP
