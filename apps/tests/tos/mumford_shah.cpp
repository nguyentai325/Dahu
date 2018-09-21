#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/morpho/tos/tos.hpp>
#include <mln/core/grays.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/io/imread.hpp>

#ifndef MLN_IMG_PATH
# error "MLN_IMG_PATH must be defined."
#endif


#define BOOST_TEST_MODULE TOS
#include <boost/test/unit_test.hpp>

#include <apps/tos/mumford_shah.hpp>


BOOST_AUTO_TEST_CASE(mumford_shah_gray)
{
  using namespace mln;

  image2d<uint8> ima;

  io::imread(MLN_IMG_PATH "/squares.pgm", ima);

   typedef UInt<9> V;
   image2d<V> K;
   std::vector<unsigned> S;
   image2d<unsigned> parent;

   image2d<V> f = transform(ima, [](uint8 x) -> V { return x*2; });
   std::tie(K, parent, S) = morpho::ToS(f, c4);

   image2d<float> curv = curvature_on_edge(ima);
   image2d< internal::energy_t<uint8> > attributes;
   image2d<float> energy = compute_energy(ima, K, parent, S, attributes);

   for (unsigned x: S)
     {
       if (K[parent[x]] != K[x] or x == parent[x])
	 std::cout << x << " : " << energy[x] << std::endl
		   << attributes[x] << std::endl;
     }

   // check attributes of each shapes (computed by hand)
   internal::energy_t<uint8> shp[7];

   shp[0].m_e_length = 128 * 4;
   shp[0].m_v_n_int = (128-10) * 5 * 4 + 25 * 4;
   shp[0].m_v_n_ext = 0;
   shp[0].m_v_sum_int = 0;
   shp[0].m_v_sum_int_sqr = 0;
   shp[0].m_v_sum_ext = 0;
   shp[0].m_v_sum_ext_sqr = 0;

   {
     int k = 1, w = 99, h = 97, glint = 130, glext = 0;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = 1771 * 130 + 89 * 252;
     shp[k].m_v_sum_int_sqr = 1771 * 130 * 130 + 89 * 252 * 252;
     shp[k].m_v_sum_ext = shp[k].m_v_n_ext * glext;
     shp[k].m_v_sum_ext_sqr = shp[k].m_v_sum_ext * glext;
   }

   {
     int k = 2, w = 89, h = 78, glint = 252, glext = 130;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = shp[k].m_v_n_int * glint;
     shp[k].m_v_sum_int_sqr = shp[k].m_v_sum_int * glint;
     shp[k].m_v_sum_ext = 1671 * 130;
     shp[k].m_v_sum_ext_sqr = 1671 * 130 * 130;
   }

   {
     int k = 3, w = 25, h = 33, glint = 195, glext = 252;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = shp[k].m_v_n_int * glint;
     shp[k].m_v_sum_int_sqr = shp[k].m_v_sum_int * glint;
     shp[k].m_v_sum_ext = shp[k].m_v_n_ext * glext;
     shp[k].m_v_sum_ext_sqr = shp[k].m_v_sum_ext * glext;
   }

   {
     int k = 4, w = 21, h = 57, glext = 252;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = 172 * 251 + 508 * 154;
     shp[k].m_v_sum_int_sqr = 172 * 251 * 251 + 508 * 154 * 154;
     shp[k].m_v_sum_ext = shp[k].m_v_n_ext * glext;
     shp[k].m_v_sum_ext_sqr = shp[k].m_v_sum_ext * glext;
   }

   {
     int k = 5, w = 17, h = 25, glint = 251;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = shp[k].m_v_n_int * glint;
     shp[k].m_v_sum_int_sqr = shp[k].m_v_sum_int * glint;
     shp[k].m_v_sum_ext = 268 * 154 + 252 * 252;
     shp[k].m_v_sum_ext_sqr = 268 * 154 * 154 + 252 * 252 * 252;
   }

   {
     int k = 6, w = 10, h = 11, glint = 12, glext=154;
     shp[k].m_e_length = w * 2 + h * 2;
     shp[k].m_v_n_int = (w-10) * 2 * 5 + (h-10) * 2 * 5 + 25 * 4;
     shp[k].m_v_n_ext = w * 2 * 5 + h * 2 * 5 + 25 * 4;
     shp[k].m_v_sum_int = shp[k].m_v_n_int * glint;
     shp[k].m_v_sum_int_sqr = shp[k].m_v_sum_int * glint;
     shp[k].m_v_sum_ext = shp[k].m_v_n_ext * glext;
     shp[k].m_v_sum_ext_sqr = shp[k].m_v_sum_ext * glext;
   }


   std::cout << "==========================" << std::endl;
   for (auto& a : shp)
     std::cout << a << std::endl;
}
