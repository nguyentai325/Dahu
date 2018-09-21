# include <mln/core/image/image2d.hpp>

# include <mln/io/imread.hpp>
# include <mln/io/imsave.hpp>

# include <mln/core/algorithm/paste.hpp>

#include <apps/tos/Kinterpolate.hpp>
#include <apps/tos/colorToSGrad.hpp>
#include <apps/tos/mumford_shah.hpp>
#include <apps/tos/addborder.hpp>
#include <apps/tos/set_mean_on_nodes.hpp>
#include <apps/tos/objdetection.hpp>
#include <boost/format.hpp>

namespace mln
{
  /*
  inline
  unsigned
  zfindroot(image2d<unsigned>& par, unsigned x)
  {
    if (par[x] != x)
      par[x] = zfindroot(par, par[x]);
    return par[x];
  }


  template <typename V>
  void
  getObjects(const image2d<float>& energy,
	     image2d<V>& K,
	     image2d<unsigned>& parent,
	     const std::vector<unsigned>& S,
	     unsigned seuil)
  {
    static const unsigned NullNode = -1;

    struct child_t {
      unsigned first_child = -1;
      unsigned next_sibling = -1;
    };

    trace::entering("Meaningfull line detection");

    // Copy energy of nodes and set child relation
    std::vector< std::pair<unsigned, float> > E;
    image2d<child_t> crel;
    resize(crel, K);

    for (unsigned x: S)
      if (K[parent[x]] != K[x]) {
	E.push_back(std::make_pair(x, energy[x]));
	unsigned y = crel[parent[x]].first_child;
	if (y == NullNode)
	  crel[parent[x]].first_child = x;
	else  {
	  while (crel[y].next_sibling != NullNode)
	    y = crel[y].next_sibling;
	  crel[y].next_sibling = x;
	}
      }

    // Sort
    std::sort(E.begin(), E.end(), [](const std::pair<unsigned, float>& a, const std::pair<unsigned, float>& b) {
	return a.second < b.second; });


    // Union-find
    unsigned n = E.size();
    image2d<unsigned> par;
    image2d<unsigned> area;
    resize(par, K);
    resize(area, K);

    unsigned num_node_deleted = 0;
    unsigned num_node_kept = 0;

    for (unsigned i = 0; i < n; ++i)
      {
	unsigned x = E[i].first;
	par[x] = x;
	area[x] = 1;

	// Visist Nbh
	// Parent First
	{
	  unsigned y = zfindroot(par, parent[x]);
	  if (y != x) {
	    par[y] = x;
	    area[x] += area[y];
	  }
	}
	// Children next
	{
	  unsigned y = crel[x].first_child;
	  if (y != NullNode)
	    do {
	      unsigned r = zfindroot(par, y);
	      if (r != x) {
		par[r] = x;
		area[x] += area[r];
	      }
	      y = crel[y].next_sibling;
	    } while (y != NullNode);
	}
	if (area[x] < seuil)
	  num_node_deleted++;
	else
	  num_node_kept++;
      }

    // simplify, remove non significant level lines
    // AND RECANONIZE !
    {
      for (unsigned x: S) {
	if (area[x] < seuil) {
	  K[x] = K[parent[x]];
	}
	unsigned q = parent[x];
	if (K[q] == K[parent[q]])
	  parent[x] = parent[q];
      }
    }

    std::cout << "Level line removed: " << num_node_deleted << std::endl
	      << "Level line kept: " << num_node_kept << std::endl;

    trace::exiting();
  }
  */

} // end of mln


namespace mln
{

  template <typename V>
  image2d<unsigned>
  make_unique_id(const image2d<V>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    image2d<unsigned> out;
    resize(out, K);
    for (unsigned x: S)
      if (K[x] == K[parent[x]])
	out[x] = parent[x];
      else
	out[x] = x;
    return out;
  }

}
void usage(char** argv)
{
  std::cout << "Usage: " << argv[0] << " input(color) output[wo ext] lambda [lambda_2...]" << std::endl
	    << "Perform a simplification of the ToS by removing non significant"
	    << "level lines. A closing in shape space in performed." << std::endl
	    << "Param:" << std::endl
	    << "lambda: minimal grain size (in shape space)" << std::endl
	    << std::endl;
  std::terminate();
}


int main(int argc, char** argv)
{
  static const bool use_tos = true;

  if (argc < 4)
    usage(argv);

  using namespace mln;
  image2d<rgb8> ima;
  io::imread(argv[1], ima);


  //auto f = interpolate_median(ima, UInt<9> ());


  image2d<rgb8> ima_ = addborder(ima, lexicographicalorder_less<rgb8>() );
  image2d<rgb8> f = interpolate_k1(ima_);

  image2d<unsigned> K;
  std::vector<unsigned> S;
  image2d<unsigned> parent;

  image2d<rgb8> out, final;
  resize(out, K);
  resize(final, ima);


  // case 1 with tos
  if (!use_tos)
    {
      colorToSGrad_with_mintree(ima, K, parent, S);

      image2d<float> energy = compute_energy(f, K, parent, S);
      image2d<unsigned> ids = make_unique_id(K, parent, S);

      for (int i = 3; i < argc; ++i)
	{
	  int lambda = std::atoi(argv[i]);
	  auto cIds = clone(ids);
	  auto cParent = clone(parent);
	  auto cEnergy = clone(energy);

	  getObjects(cEnergy, cIds, cParent, S, lambda);
	  out = set_mean_on_node(f, cIds, S, cParent, K1::is_face_2);

	  {
	    box2d d = out.domain();
	    sbox2d sub_domain { {2,2}, {d.pmax[0]-2, d.pmax[1]-2}, {2,2} };
	    copy(out | sub_domain, final);
	  }

	  if (argc == 4)
	    io::imsave(final, (std::string(argv[2]) + ".tiff").c_str());
	  else
	    io::imsave(final, (boost::format("%s-%06i.tiff") % argv[2] % lambda).str().c_str());
	}

    }
  else
    {
      colorToSGrad(ima, K, parent, S);
      //colorToSGrad_f2(ima, K, parent, S);

      image2d<float> energy = compute_energy(f, K, parent, S);
      //image2d<float> energy = compute_energy(ima, K, parent, S);
      image2d<unsigned> ids = make_unique_id(K, parent, S);

      for (int i = 3; i < argc; ++i)
	{
	  int lambda = std::atoi(argv[i]);
	  auto cIds = clone(ids);
	  auto cParent = clone(parent);
	  auto cEnergy = clone(energy);

	  getObjects(cEnergy, cIds, cParent, S, lambda);
	  out = set_mean_on_node(interpolate_k1(f), cIds, S, cParent, K1::is_face_2);
	  //out = set_mean_on_node(f, cIds, S, cParent, K1::is_face_2);

	  {
	    box2d d = out.domain();
	    //sbox2d sub_domain { {2,2}, {d.pmax[0]-2, d.pmax[1]-2}, {2,2} };
	    sbox2d sub_domain { {4,4}, {d.pmax[0]-4, d.pmax[1]-4}, {4,4} };
	    copy(out | sub_domain, final);
	  }

	  if (argc == 4)
	    io::imsave(final, (std::string(argv[2]) + ".tiff").c_str());
	  else
	    io::imsave(final, (boost::format("%s-%06i.tiff") % argv[2] % lambda).str().c_str());
	}
    }
}
