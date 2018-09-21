#ifndef PIX_RANGE_ITERATOR_PROXY_HPP
# define PIX_RANGE_ITERATOR_PROXY_HPP

# include <boost/range/iterator.hpp>

namespace mln
{

  namespace internal
  {

    template <typename PixelRange>
    struct pix_range_iterator_proxy
    {
      typedef typename boost::range_iterator<PixelRange>::type iterator;

      pix_range_iterator_proxy(const PixelRange& r);

      void init();
      void next();
      bool is_valid() const;

      // As a proxy to pixel
      typedef typename std::iterator_traits<iterator>::value_type pixel_t;
      typedef typename std::iterator_traits<iterator>::reference  ref_pixel_t;
      typedef typename pixel_t::site_type       site_type;
      typedef typename pixel_t::point_type      point_type;
      typedef typename pixel_t::value_type      value_type;
      typedef typename pixel_t::reference       reference;
      typedef typename pixel_t::image_type      image_type;

      static_assert(std::is_reference<ref_pixel_t>::value,
                    "It is required that pixel iterator returns real reference to a pixel object.");


      operator ref_pixel_t () const      { return *cur_; }
      reference   val() const            { return cur_->val(); }
      site_type   site() const           { return cur_->site(); }
      point_type  point() const          { return cur_->point(); }
      image_type& image() const          { return cur_->image(); }


      /// To use in a forall loop
      /// \{
      std::false_type set_dejavu_(bool v);
      bool get_dejavu_() const;
      /// \}

    private:
      pix_range_iterator_proxy(const pix_range_iterator_proxy& other);

    private:
      PixelRange r_;
      iterator cur_;
      iterator end_;
      bool deja_vu_;
    };



    /***********/
    /* Implem  */
    /***********/


    template <typename PixelRange>
    pix_range_iterator_proxy<PixelRange>::pix_range_iterator_proxy(const PixelRange& r) :
      r_ (r), cur_ (r_.begin()), end_(r_.end())
    {
    }

    template <typename PixelRange>
    void
    pix_range_iterator_proxy<PixelRange>::init()
    {
      cur_ = r_.begin();
      end_ = r_.end();
      deja_vu_ = false;
    }

    template <typename PixelRange>
    void
    pix_range_iterator_proxy<PixelRange>::next()
    {
      if (!deja_vu_)
	++cur_;
      deja_vu_ = true;
    }

    template <typename PixelRange>
    bool
    pix_range_iterator_proxy<PixelRange>::is_valid() const
    {
      return cur_ != end_;
    }


    template <typename PixelRange>
    std::false_type
    pix_range_iterator_proxy<PixelRange>::set_dejavu_(bool v)
    {
      deja_vu_ = v;
      return std::false_type ();
    }

    template <typename PixelRange>
    bool
    pix_range_iterator_proxy<PixelRange>::get_dejavu_() const
    {
      return deja_vu_;
    }
  }

}

#endif // ! PIX_RANGE_ITERATOR_PROXY_HPP
