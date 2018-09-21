#ifndef MLN_CORE_SLIDING_OBJECT_HPP
# define MLN_CORE_SLIDING_OBJECT_HPP
namespace mln {

  namespace internal
  {
    template <typename Image, typename image_tag>
    struct sliding_object_base

  }

  template <typename Image>
  struct sliding_object
  {

  private:
    Image& ima;
    
  }

  template <typename Image, typename SiteSet>
  struct meta_sliding_object
  {
    meta_sliding_object(const Image& ima, const SiteSet& dpoints)
    {
      ima
    }

  };


} // end of namespace mln



#endif //!MLN_CORE_SLIDING_OBJECT_HPP
