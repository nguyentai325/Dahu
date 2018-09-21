#ifndef MLN_CORE_VEC_HPP
# define MLN_CORE_VEC_HPP

# include <mln/core/value/value_traits.hpp>
# include <mln/core/vec_base.hpp>
# include <mln/core/vec/compare.hpp>
# include <mln/core/vec/vec_io.hpp>
# include <mln/core/vec/vec_ops.hpp>
# include <mln/core/vec/vec_math_ops.hpp>


namespace mln
{

  template <typename T, unsigned dim>
  using vec = internal::vec_base<T, dim, generic_vector_tag>;

  typedef vec<unsigned char, 1>	vec1ub;
  typedef vec<unsigned char, 2>	vec2ub;
  typedef vec<unsigned char, 3>	vec3ub;
  typedef vec<char, 1>		vec1b;
  typedef vec<char, 2>		vec2b;
  typedef vec<char, 3>		vec3b;
  typedef vec<unsigned short, 1> vec1us;
  typedef vec<unsigned short, 2> vec2us;
  typedef vec<unsigned short, 3> vec3us;
  typedef vec<short, 1>		vec1s;
  typedef vec<short, 2>		vec2s;
  typedef vec<short, 3>		vec3s;
  typedef vec<unsigned, 1>	vec1u;
  typedef vec<unsigned, 2>	vec2u;
  typedef vec<unsigned, 3>	vec3u;
  typedef vec<int, 1>		vec1i;
  typedef vec<int, 2>		vec2i;
  typedef vec<int, 3>		vec3i;
  typedef vec<float, 1>		vec1f;
  typedef vec<float, 2>		vec2f;
  typedef vec<float, 3>		vec3f;

}

#endif // ! MLN_CORE_VEC_HPP
