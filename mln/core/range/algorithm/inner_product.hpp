#ifndef MLN_CORE_RANGE_ALGORITHM_INNER_PRODUCT_HPP
# define MLN_CORE_RANGE_ALGORITHM_INNER_PRODUCT_HPP

# include <mln/core/range/range.hpp>


namespace mln
{
  namespace range
  {

    template <typename Range1, typename Range2, typename T>
    T inner_product(const Range1& r1, const Range2& r2, T value)
    {
      mln_iter(vin1, r1);
      mln_iter(vin2, r2);
      mln_forall(vin1, vin2)
        value = value + *vin1 * *vin2;
      return value;
    }

    template <typename Range1, typename Range2, typename T, class BinaryOperation1, class BinaryOperation2>
    T inner_product(const Range1& r1, const Range2& r2, T value, BinaryOperation1 op1, BinaryOperation2 op2)
    {
      mln_iter(vin1, r1);
      mln_iter(vin2, r2);
      mln_forall(vin1, vin2)
        value = op2(value, op1(*vin1,*vin2));
      return value;
    }

  } // end of namespace mln::range
} // end of namespace mln

#endif //!MLN_CORE_RANGE_ALGORITHM_INNER_PRODUCT_HPP
