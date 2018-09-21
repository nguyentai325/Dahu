#ifndef MLN_DATA_STRETCH_HPP
# define MLN_DATA_STRETCH_HPP

# include <mln/core/image/image.hpp>
# include <mln/core/algorithm/transform.hpp>
# include <mln/core/algorithm/accumulate.hpp>
# include <mln/accu/accumulators/minmax.hpp>

namespace mln
{

  namespace data
  {

    /// \brief Transform image values to fit the complete dynamic of \p V
    ///
    /// If V is integral, values are streched to \f$[V_{min}, V_{max}]\f$, otherwise,
    /// if V is floating point, values are streched to \f$[0,1]\f. The method only support
    /// V and not vectorial type.
    template <class V, class I>
    mln_ch_value(I, V) stretch(const Image<I>& f);


    template <class I, class J>
    J& stretch_to(const Image<I>& f, Image<J>& out);

    template <class I, class J>
    J&& stretch_to(const Image<I>& f, Image<J>&& out);

    /************************/
    /*** Implementation   ***/
    /************************/

    template <class I, class J>
    J& stretch_to(const Image<I>& f, Image<J>& out)
    {
      typedef mln_value(J) V;

      static_assert(std::is_convertible<mln_value(I), V>::value,
                    "The image value type must be convertible to V");
      static_assert(std::is_arithmetic<V>::value,
                    "V must be arithmetic");

      mln_entering("mln::data:stretch");

      double m, M;
      std::tie(m, M) = accumulate(f, accu::accumulators::minmax<mln_value(I)> ());

      double x = (not std::is_floating_point<V>::value) ? (double)value_traits<V>::min() : 0.0d;
      double y = (not std::is_floating_point<V>::value) ? (double)value_traits<V>::max() : 1.0d;

      double r = (y - x) / (M - m);
      transform(f, [m,x,r](mln_value(I) v) -> V { return x + (v - m) * r; }, out);
      mln_exiting();
      return exact(out);
    }

    template <class I, class J>
    J&& stretch_to(const Image<I>& f, Image<J>&& out)
    {
      stretch_to(f, out);
      return move_exact(out);
    }


    template <class V, class I>
    mln_ch_value(I, V) stretch(const Image<I>& f)
    {
      static_assert(std::is_convertible<mln_value(I), V>::value,
                    "The image value type must be convertible to V");
      static_assert(std::is_arithmetic<V>::value,
                    "V must be arithmetic");

      mln_entering("mln::data:stretch");

      double m, M;
      std::tie(m, M) = accumulate(f, accu::accumulators::minmax<mln_value(I)> ());

      double x = (not std::is_floating_point<V>::value) ? (double)value_traits<V>::min() : 0.0d;
      double y = (not std::is_floating_point<V>::value) ? (double)value_traits<V>::max() : 1.0d;

      double r = (y - x) / (M - m);
      mln_ch_value(I, V) out = transform(f, [m,x,r](mln_value(I) v) -> V { return x + (v - m) * r; });
      mln_exiting();

      return out;
    }


  }

}

#endif // ! MLN_DATA_STRETCH_HPP
