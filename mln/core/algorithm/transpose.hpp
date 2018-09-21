#ifndef MLN_CORE_ALGORITHM_TRANSPOSE_HPP
# define MLN_CORE_ALGORITHM_TRANSPOSE_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/trace.hpp>

namespace mln
{

  template <class I, class J>
  void
  transpose(const Image<I>& ima,
            Image<J>& out);

  template <class I>
  image2d<mln_value(I)>
  transpose(const Image<I>& ima);

  /******************************************/
  /****          Implementation          ****/
  /******************************************/

  template <class I, class J>
  void
  transpose(const Image<I>& ima_,
            Image<J>& out_)
  {
    mln_entering("mln::transpose");

    const I& ima = exact(ima_);
    J& out = exact(out_);

    mln_foreach(point2d p, ima.domain())
      out(point2d{p[1],p[0]}) = ima(p);

    mln_exiting();
  }

  namespace internal
  {

    template <class I>
    image2d<mln_value(I)>
    transpose_init_out(const I& ima, box2d d, extension::extension_tag)
    {
      (void) ima;
      return image2d<mln_value(I)>(d);
    }

    template <class I>
    image2d<mln_value(I)>
    transpose_init_out(const I& ima, box2d d, extension::border_extension_tag)
    {
      return image2d<mln_value(I)>(d, ima.border());
    }
  }


  template <class I>
  image2d<mln_value(I)>
  transpose(const Image<I>& ima_)
  {
    const I& ima = exact(ima_);

    box2d dom = ima.domain();
    box2d dom2 = {{dom.pmin[1], dom.pmin[0]},
                  {dom.pmax[1], dom.pmax[0]}};

    image2d<mln_value(I)> out = internal::transpose_init_out(ima, dom2, typename image_traits<I>::extension ());
    transpose(ima, out);
    return out;
  }


} // end of namespace mln

#endif //!MLN_CORE_ALGORITHM_TRANSPOSE_HPP
