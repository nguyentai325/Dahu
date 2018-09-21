#ifndef MLN_CORE_ITERATOR_SLIDING_WIN_PIXTER_HPP
# define MLN_CORE_ITERATOR_SLIDING_WIN_PIXTER_HPP

# include <mln/core/iterator/iterator_base.hpp>
# include <mln/core/std/array.hpp>
# include <mln/core/wrt_offset.hpp>
# include <mln/core/range/size.hpp>
# include <vector>

namespace mln
{


  template <typename SiteSet, typename PixelOrPixelIterator>
  struct sliding_win_pixter;


  namespace internal
  {
    template <typename P, typename V, typename I, typename _category = typename image_traits<I>::category>
    struct sliding_win_pixel;


    template <typename SiteSet, typename PixelOrPixelIterator, typename _category, typename E>
    struct sliding_win_pixter_base;

    template <typename SiteSet, typename PixelOrPixelIterator, typename E>
    struct sliding_win_pixter_base< SiteSet, PixelOrPixelIterator, raw_image_tag, E>;

    template <typename Point, size_t N, typename PixelOrPixelIterator, typename E>
    struct sliding_win_pixter_base< mln::array<Point, N>, PixelOrPixelIterator, raw_image_tag, E>;

    template <typename PixelOrPixelIterator, bool is_iterator = is_a<PixelOrPixelIterator, Iterator>::value >
    struct sliding_win_pixter_dispatch;
  }

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  namespace internal
  {

    template <typename PixelOrPixelIterator>
    struct sliding_win_pixter_dispatch<PixelOrPixelIterator, true>
    {
      typedef typename PixelOrPixelIterator::value_type Pixel;
      static typename PixelOrPixelIterator::reference getpixel(const PixelOrPixelIterator* x) { return x->dereference(); } ;
    };

    template <typename PixelOrPixelIterator>
    struct sliding_win_pixter_dispatch<PixelOrPixelIterator, false>
    {
      typedef PixelOrPixelIterator Pixel;
      static const Pixel& getpixel(const PixelOrPixelIterator* x) { return *x; } ;
    };


  }

  template <typename SiteSet, typename PixelOrPixelIterator>
  struct sliding_win_pixter :
    internal::sliding_win_pixter_base< SiteSet, PixelOrPixelIterator,
				       typename image_traits<typename internal::sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::image_type>::category,
				       sliding_win_pixter<SiteSet, PixelOrPixelIterator> >
  {
  private:
    typedef internal::sliding_win_pixter_base< SiteSet, PixelOrPixelIterator,
					       typename image_traits<typename internal::sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::image_type>::category,
                                               sliding_win_pixter<SiteSet, PixelOrPixelIterator> > base;

  public:
    sliding_win_pixter() = default;
    sliding_win_pixter(const SiteSet& domain, const PixelOrPixelIterator& pix)
      : base(domain, pix)
    {
    }

  };

  /******************************************/
  /****             Gneneric             ****/
  /******************************************/
  namespace internal
  {

    template <typename P, typename V, typename I, typename _category>
    struct sliding_win_pixel
    {
      static_assert(image_traits<I>::accessible::value, "Cannot use an neighborhood/SE/win on a non-accessible image.");

      typedef P point_type;
      typedef P site_type;
      typedef V value_type;
      typedef V& reference;
      typedef I image_type;

      sliding_win_pixel() = default;
      sliding_win_pixel(const sliding_win_pixel&) = default;

      // Non-const -> const conversion
      template <typename U, typename J>
      sliding_win_pixel(const sliding_win_pixel<P, U, J>& x,
                        typename std::enable_if< std::is_convertible<typename std::remove_reference<U>::type*,
                        typename std::remove_reference<V>::type*>::value >::type* dummy = NULL)
	: m_p (x.m_p), m_ima (x.m_ima)
      {
      }


      V& val() const { return m_ima->at(m_p); }
      const P& point() const { return m_p; }
      const P& site() const { return m_p; }
      image_type& image() const { return *m_ima; }

    private:
      template <typename, typename, typename, typename>
      friend struct sliding_win_pixter_base;

      P  m_p;
      I* m_ima;
    };


    template <typename SiteSet, typename PixelOrPixelIterator, typename _category, typename E>
    struct sliding_win_pixter_base :
      iterator_base<E, const sliding_win_pixel<typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::point_type,
                                               typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::value_type,
                                               typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::image_type> >,
      protected sliding_win_pixter_dispatch<PixelOrPixelIterator>
    {
    private:
      typedef typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel Px;
      typedef sliding_win_pixel<typename Px::point_type, typename Px::value_type, typename Px::image_type> pixel_t;

    public:
      sliding_win_pixter_base() = default;

      sliding_win_pixter_base(const SiteSet& domain, const PixelOrPixelIterator& pix)
	: m_bind(&pix), m_pit(domain.iter())
      {
        m_pix.m_ima = &this->getpixel(m_bind).image();
      }

      void init()
      {
        m_pit.init();
        m_pix.m_p = this->getpixel(m_bind).point() + *m_pit;
      }

      void next()
      {
        m_pit.next();
        m_pix.m_p = this->getpixel(m_bind).point() + *m_pit;
      }

      bool finished() const
      {
        return m_pit.finished();
      }

      const pixel_t& dereference() const
      {
        return m_pix;
      }


    private:
      const PixelOrPixelIterator*                       m_bind;
      typename range_iterator<const SiteSet>::type      m_pit;

      pixel_t m_pix;
    };



  }



  /******************************************/
  /****          Specialization          ****/
  /******************************************/

  namespace internal
  {

    template <typename P, typename V, typename I>
    struct sliding_win_pixel<P, V, I, raw_image_tag>
    {
      typedef P point_type;
      typedef P site_type;
      typedef V value_type;
      typedef V& reference;
      typedef I image_type;

      sliding_win_pixel() = default;
      sliding_win_pixel(const sliding_win_pixel&) = default;

      // Non-const -> const conversion
      template <typename U, typename J>
      sliding_win_pixel(const sliding_win_pixel<P, U, J>& x,
                        typename std::enable_if< std::is_convertible<typename std::remove_reference<U>::type*,
                        typename std::remove_reference<V>::type*>::value >::type* dummy = NULL)
	: v_ (x.v_), p_ (x.p_), image_ (x.image_)
      {
      }


      V& val() const { return *v_; }
      const P& point() const { return p_; }
      const P& site() const { return p_; }
      image_type& image() const { return *image_; }
      unsigned index() const { return idx_; }

    private:
      template <typename, typename, typename, typename>
      friend struct sliding_win_pixter_base;

      V* v_;
      P  p_;
      unsigned idx_;
      I* image_;
    };

    template <typename SiteSet, typename PixelOrPixelIterator, typename E>
    struct sliding_win_pixter_base< SiteSet, PixelOrPixelIterator, raw_image_tag, E>
      : iterator_base<E, const sliding_win_pixel<typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::point_type,
                                                 typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::value_type,
                                                 typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::image_type> >,
        protected sliding_win_pixter_dispatch<PixelOrPixelIterator>
    {
    private:
      typedef typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel Pixel;
      typedef typename Pixel::point_type P;
      typedef typename Pixel::image_type I;
      typedef typename Pixel::value_type V;
      typedef typename I::difference_type difference_type;
      typedef sliding_win_pixel<P,V,I> pixel_t;

    public:
      sliding_win_pixter_base() = default;

      sliding_win_pixter_base(const SiteSet& domain, const PixelOrPixelIterator& pix)
	: bind_ (&pix), domain_(&domain)
      {
	unsigned sz = rng::size(domain);
	offset_.resize(sz);
	indexes_.resize(sz);
	wrt_offset( this->getpixel(bind_).image(), domain, offset_.begin());
	wrt_delta_index( this->getpixel(bind_).image(), domain, indexes_.begin());
	pit_ = domain.iter();
      }

      void init() {
	i_ = 0;
	p_ = this->getpixel(bind_).point();
	idx_ = this->getpixel(bind_).index();
	ptr_ = (char*)(& this->getpixel(bind_).val());
	pit_.init();
	mypix_.p_ = p_ + *pit_;
	mypix_.v_ = (V*)(ptr_ + offset_[0]);
	mypix_.idx_ = idx_ + indexes_[0];
      }

      void next() {
	++i_;
	pit_.next();
	mypix_.p_ = p_ + *pit_;
	mypix_.v_ = (V*)(ptr_ + offset_[i_]);
	mypix_.idx_ = idx_ + indexes_[i_];
      }

      bool finished() const {
	return i_ == offset_.size();
      }

      const pixel_t&
      dereference() const {
	return mypix_;
      }

    private:
      const PixelOrPixelIterator* bind_;
      const SiteSet* domain_;
      std::vector<difference_type> offset_;
      std::vector<difference_type> indexes_;
      typename SiteSet::iterator pit_;
      pixel_t mypix_;

      unsigned i_;
      unsigned idx_;
      char* ptr_;
      P p_;

    };



    ////////// Specialization for SiteSet = mln::array   ///////////////

    template <typename Point, size_t N, typename PixelOrPixelIterator, typename E>
    struct sliding_win_pixter_base< mln::array<Point, N>, PixelOrPixelIterator, raw_image_tag, E>
      : iterator_base<E, const sliding_win_pixel<typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::point_type,
                                                 typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::value_type,
                                                 typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel::image_type> >,
	protected sliding_win_pixter_dispatch<PixelOrPixelIterator>
    {
    private:
      typedef typename sliding_win_pixter_dispatch<PixelOrPixelIterator>::Pixel Pixel;
      typedef typename Pixel::point_type P;
      typedef typename Pixel::image_type I;
      typedef typename Pixel::value_type V;
      typedef typename I::difference_type difference_type;
      typedef sliding_win_pixel<P,V,I> pixel_t;

    public:
      sliding_win_pixter_base() = default;

      sliding_win_pixter_base(const mln::array<Point, N>& domain, const PixelOrPixelIterator& pix)
	: bind_ (&pix), domain_(&domain)
      {
	wrt_offset( this->getpixel(bind_).image(), domain, offset_);
	wrt_delta_index( this->getpixel(bind_).image(), domain, indexes_);
      }

      void init() {
	i_ = 0;
	idx_ = this->getpixel(bind_).index();
	p_ = this->getpixel(bind_).point();
	ptr_ = (char*)(& this->getpixel(bind_).val());
	mypix_.p_ = p_ + (*domain_)[0];
	mypix_.v_ = (V*)(ptr_ + offset_[0]);
	mypix_.idx_ = idx_ + indexes_[0];
      }

      void next() {
	++i_;
	mypix_.p_ = p_ + (*domain_)[i_];
	mypix_.v_ = (V*)(ptr_ + offset_[i_]);
	mypix_.idx_ = idx_ + indexes_[i_];
      }

      bool finished() const {
	return i_ == N;
      }

      const pixel_t&
      dereference() const {
	return mypix_;
      }

    private:
      const PixelOrPixelIterator* bind_;
      const mln::array<Point, N>* domain_;
      mln::array<difference_type, N> offset_;
      mln::array<difference_type, N> indexes_;
      pixel_t mypix_;

      int i_;
      char* ptr_;
      Point p_;
      unsigned idx_;

    };



  }

} // end of namespace mln


#endif //!MLN_CORE_ITERATOR_SLIDING_WIN_PIXTER_HPP
