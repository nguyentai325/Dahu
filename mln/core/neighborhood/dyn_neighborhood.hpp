#ifndef MLN_CORE_NEIGHBORHOOD_DYN_NEIGHBORHOOD_HPP
# define MLN_CORE_NEIGHBORHOOD_DYN_NEIGHBORHOOD_HPP

# include <mln/core/neighborhood/neighborhood_base.hpp>
# include <mln/core/range/range.hpp>
# include <mln/core/range/iterator_range.hpp>
# include <mln/core/neighborhood/sliding_pixter.hpp>
# include <mln/core/neighborhood/sliding_piter.hpp>


namespace mln
{

  template <class SiteSet, class category, class E>
  struct dyn_neighborhood_base;

  template <class SiteSet, class category>
  struct dyn_neighborhood;

  /******************************************/
  /****          Implementation          ****/
  /******************************************/


  // template <class SiteSet, class category, class E>
  // struct dyn_neighborhood_base
  // {
  //   dyn_neighborhood_base(const SiteSet& x)
  //   {
  //   }
  // };


  template <class SiteSet, class E>
  struct dyn_neighborhood_base<SiteSet, dynamic_neighborhood_tag, E>
    : neighborhood_base<E, dynamic_neighborhood_tag>
  {
  public:
    typedef typename range_value<SiteSet>::type point_type;
    typedef typename range_value<SiteSet>::type site_type;

    dyn_neighborhood_base(const SiteSet& dpts)
    : dpoints(dpts)
    {
      m_dpoints_acc.resize(rng::size(dpoints));
      mln_iter(p, dpoints);
      p.init();
      m_dpoints_acc[0] = *p;
      for (int i = 1; p.finished(); ++i, p.next())
	m_dpoints_acc[i] = *p - m_dpoints_acc[i-1];
    }

    ~dyn_neighborhood_base() = default;

    iterator_range< sliding_piter<point_type, SiteSet> >
    __process_point(const point_type& p) const
    {
      return make_iterator_range( sliding_piter<point_type, SiteSet>(p, this->dpoints) );
    }

    template <typename P>
    iterator_range< sliding_piter<const P*, SiteSet> >
    __bind_point(P& p) const
    {
      return make_iterator_range( sliding_piter<const P*, SiteSet>(&p, this->dpoints) );
    }

    template <typename PointIterator>
    iterator_range< sliding_piter<PointIterator, SiteSet> >
    __bind_point_iterator(const PointIterator& p) const
    {
      return make_iterator_range( sliding_piter<PointIterator, SiteSet>(p, this->dpoints) );
    }

    template <typename Px>
    iterator_range< sliding_pixter<const Px*, SiteSet> >
    __bind_pixel(Px& px) const
    {
      return make_iterator_range( sliding_pixter<const Px*, SiteSet>(&px, this->dpoints) );
    }

    template <typename Px>
    iterator_range< sliding_pixter<Px, SiteSet> >
    __bind_pixel_iterator(const Px& px) const
    {
      return make_iterator_range( sliding_pixter<Px, SiteSet>(px, this->dpoints) );
    }

    const SiteSet dpoints;

  private:
    std::vector<point_type> m_dpoints_acc;
  };


  template <class SiteSet, class category>
  struct dyn_neighborhood
    : dyn_neighborhood_base<SiteSet,
			    category,
			    dyn_neighborhood<SiteSet, category> >
  {
    dyn_neighborhood(const SiteSet& s)
    : dyn_neighborhood_base<SiteSet,
			    category,
			    dyn_neighborhood<SiteSet, category> >(s)
    {
    }
  };


  /**************************************/
  /***   Specialization of static  NBH **/
  /**************************************/

  template <class SiteSet, class E>
  struct dyn_neighborhood_base<SiteSet, constant_neighborhood_tag, E>
    : neighborhood_base<E, constant_neighborhood_tag>
  {
  private:
    enum { N = rng::size(*(SiteSet*)0) };

  public:
    typedef typename range_value<SiteSet>::type point_type;
    typedef typename range_value<SiteSet>::type site_type;


    dyn_neighborhood_base()
    {
      if (not __init__)
	{
	  const SiteSet& dpoints = exact(this)->dpoints;
	  m_dpoints_acc[0] = dpoints[0];
	  for (int i = 1; i < N; ++i)
	    m_dpoints_acc[i] = dpoints[i] - m_dpoints_acc[i-1];
	  __init__ = true;
	}
    }

    ~dyn_neighborhood_base() = default;

    iterator_range< sliding_piter<point_type, SiteSet> >
    __process_point(const point_type& p) const
    {
      return make_iterator_range( sliding_piter<point_type, SiteSet>(p, exact(this)->dpoints) );
    }

    template <typename P>
    iterator_range< sliding_piter<const P*, SiteSet> >
    __bind_point(P& p) const
    {
      return make_iterator_range( sliding_piter<const P*, SiteSet>(&p, exact(this)->dpoints) );
    }

    template <typename PointIterator>
    iterator_range< sliding_piter<PointIterator, SiteSet> >
    __bind_point_iterator(const PointIterator& p) const
    {
      return make_iterator_range( sliding_piter<PointIterator, SiteSet>(p, exact(this)->dpoints) );
    }

    template <typename Px>
    iterator_range< sliding_pixter<const Px*, SiteSet> >
    __bind_pixel(Px& px) const
    {
      return make_iterator_range( sliding_pixter<const Px*, SiteSet>(&px, exact(this)->dpoints) );
    }

    template <typename Px>
    iterator_range< sliding_pixter<Px, SiteSet> >
    __bind_pixel_iterator(const Px& px) const
    {
      return make_iterator_range( sliding_pixter<Px, SiteSet>(px, exact(this)->dpoints) );
    }


  private:
    static bool __init__;
    static std::array<point_type, N> m_dpoints_acc;
  };

  template <class SiteSet, class E>
  bool dyn_neighborhood_base<SiteSet, constant_neighborhood_tag, E>::__init__ = false;

  template <class SiteSet, class E>
  std::array<typename range_value<SiteSet>::type,
	     dyn_neighborhood_base<SiteSet, constant_neighborhood_tag, E>::N>
  dyn_neighborhood_base<SiteSet, constant_neighborhood_tag, E>::m_dpoints_acc;


  template <class SiteSet>
  struct dyn_neighborhood<SiteSet, constant_neighborhood_tag>
    : dyn_neighborhood_base<SiteSet,
			    constant_neighborhood_tag,
			    dyn_neighborhood<SiteSet, constant_neighborhood_tag> >
  {
    dyn_neighborhood(const SiteSet& s)
    : dpoints(s)
    {
    }

    SiteSet dpoints;
  };


}

#endif // ! MLN_CORE_NEIGHBORHOOD_DYN_NEIGHBORHOOD_HPP
