#ifndef MLN_CORE_RANGE_ADAPTOR_TRANSFORMED_HPP
# define MLN_CORE_RANGE_ADAPTOR_TRANSFORMED_HPP

# include <mln/core/range/range.hpp>
# include <mln/core/iterator/transform_iterator.hpp>

namespace mln
{

  namespace adaptors
  {

    template <typename Range, typename Function>
    struct transformed_range
    {
    private:
      typedef typename std::remove_reference<Range>::type R;

    public:
      typedef transform_iterator<typename range_iterator<R>::type, Function> iterator;
      typedef transform_iterator<typename range_const_iterator<R>::type, Function> const_iterator;

      transformed_range(Range rng, const Function& f)
        : rng_ (rng), f_ (f)
      {
      }

      // transformed_range(const transformed_rangeRange rng, const Function& f)
      // {
      // }

      iterator iter()
      {
        return iterator(rng_.iter(), f_);
      }

      const_iterator iter() const
      {
        return const_iterator(rng_.iter(), f_);
      }

    private:
      Range rng_;
      Function f_;
    };


    template <typename Range, typename Function>
    transformed_range<Range&, Function>
    transform(Range& rng, const Function& f)
    {
      return transformed_range<Range&, Function>(rng, f);
    }

  } // end of namespace mln::adaptors

} // end of namespace mln



#endif //!MLN_CORE_RANGE_ADAPTOR_TRANSFORMED_HPP
