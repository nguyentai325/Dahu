#include <mln/core/image/image2d.hpp>
#include <mln/core/neighb2d.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/extension/fill.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>
#include <boost/format.hpp>
#include "addborder.hpp"
#include "topology.hpp"

namespace mln
{


  template <typename W>
  struct myheap
  {
  private:
    static constexpr unsigned UNDEF = value_traits<unsigned>::max();

  public:
    myheap(const image2d<W>& weights)
      : m_w(weights)
    {
      resize(m_pos, weights).init(UNDEF);
      m_heap.reserve(weights.domain().size());
    }

    void push(const point2d& p)
    {
      mln_precondition(m_pos(p) == UNDEF);
      m_heap.push_back(p);
      m_pos(p) = m_heap.size() - 1;
      update(p);
    }

    point2d pop()
    {
      std::swap(m_heap.front(), m_heap.back());
      {
	point2d p = m_heap.front();
	m_pos(p) = 0;
	unsigned i = 0;
	unsigned sz = m_heap.size()-1;
	while (true)
	  {
	    unsigned minpos = i;
	    if (LCHILD(i) < sz and m_w(m_heap[minpos]) > m_w(m_heap[LCHILD(i)]))
	      minpos = LCHILD(i);
	    if (RCHILD(i) < sz and m_w(m_heap[minpos]) > m_w(m_heap[RCHILD(i)]))
	      minpos = RCHILD(i);
	    if (minpos == i)
	      break;
	    std::swap(m_heap[i], m_heap[minpos]);
	    m_pos(m_heap[i]) = i;
	    m_pos(m_heap[minpos]) = minpos;
	    i = minpos;
	  }
      }

      point2d p = m_heap.back();
      m_pos(p) = UNDEF;
      m_heap.pop_back();
      return p;
    }

    void update(const point2d& p)
    {
      mln_precondition(m_pos(p) != UNDEF);
      unsigned i = m_pos(p);
      unsigned myweight = m_w(p);
      while (i > 0 and myweight < m_w(m_heap[PAR(i)])) {
	m_heap[i] = m_heap[PAR(i)];
	m_pos(m_heap[i]) = i;
	i = PAR(i);
      }
      m_heap[i] = p;
      m_pos(p) = i;
      //std::cout << "    Position: " << i << std::endl;
      // {
      // 	std::cout << "[";
      // 	for (int i = 0; i < m_heap.size(); ++i)
      // 	  std::cout << m_heap[i] << ",";
      // 	std::cout << "]" << std::endl;;
      // }
    }

    bool empty() const { return m_heap.empty(); }

  private:
    static unsigned PAR(unsigned i) { return (i-1) / 2; }
    static unsigned LCHILD(unsigned i) { return 2*i+1; }
    static unsigned RCHILD(unsigned i) { return 2*i+2; }

    const image2d<W>&		m_w;
    std::vector<point2d>	m_heap;
    image2d<unsigned>		m_pos;
  };

  template <typename W>
  constexpr unsigned myheap<W>::UNDEF;


  unsigned
  l2norm(const rgb8& a, const rgb8& b)
  {
    unsigned v = 0;
    for (int i = 0; i < 3; ++i)
      v += std::abs( (int)a[i] - (int)b[i] );
    return v;
  }

  //
  image2d<unsigned>
  distancef(const image2d<rgb8>& f)
  {
    image2d<unsigned> grad(f.nrows()*2-1, f.ncols()*2-1);
    // Set gradient on edges
    {
      mln_foreach(point2d p, f.domain())
	{
	  point2d q = 2*p;
	  grad.at(q + point2d{0,1}) = l2norm(f(p), f.at(p + point2d{0,1}));
	  grad.at(q + point2d{1,0}) = l2norm(f(p), f.at(p + point2d{1,0}));
	}
    }

    io::imsave(transform(grad, [] (unsigned x) -> float { return x; }),
	       "gradient.tiff");
    // Dijkstra
    image2d<unsigned> mindist;

    static constexpr unsigned INFTY = value_traits<unsigned>::max();
    resize(mindist, grad).init(INFTY);

    {
      myheap<unsigned> heap(mindist);
      image2d<bool> deja_vu;

      resize(deja_vu, mindist).init(false);
      extension::fill(deja_vu, true); // Mark nodes outside domain

      point2d p{0,0};
      mindist(p) = 0;
      heap.push(p);

      std::array<point2d, 4> nbhx { {{-2,0}, {0,-2}, {0,2}, {2,0}} };

      while (! heap.empty())
	{
	  p = heap.pop();
	  //std::cout << "Popping" << p << " @ " << mindist(p) << std::endl;
	  mln_assertion(K1::is_face_2(p));
	  mln_assertion(mindist.domain().has(p));
	  deja_vu(p) = true;

	  unsigned i = 0;
	  mln_foreach(point2d e, c4(p))
	    {
	      point2d q = p + nbhx[i++];
	      if (!deja_vu.at(q))
		{
		  if (mindist(q) == INFTY) {
		    //std::cout << "  Pushing: " << q << std::endl;
		    heap.push(q);
		  }
		  unsigned old = mindist(q);
		  mindist(q) = std::min(mindist(q), mindist(p) + grad(e));
		  //std::cout << "  Updating: " << q << " from " << old << " to " << mindist(q) << std::endl;
		  heap.update(q);
		}
	    }
	}
    }

    return mindist;
  }
}

int main(int argc, char** argv)
{
  using namespace mln;

  if (argc < 3)
    {
      std::cerr << "Usage: " << argv[0] << " input.ppm output(wo extension)" << std::endl;
      std::abort();
    }

  const char* filename = argv[1];
  const char* outfile = argv[2];

  image2d<rgb8> ima;
  io::imread(filename, ima);

  image2d<rgb8> f = addborder(ima);

  image2d<unsigned> distance;
  distance = distancef(f);

  io::imsave(transform(distance, [](unsigned x) -> float { return x; }),
	     (boost::format("%s_1.tiff") % outfile).str().c_str());


  image2d<unsigned> out;
  resize(out, f);
  auto k1tok0 = sbox2d(distance.domain().pmin, distance.domain().pmax, {2,2});
  copy(distance | k1tok0, out);
  io::imsave(transform(out, [](unsigned x) -> float { return x; }),
	     (boost::format("%s_0.tiff") % outfile).str().c_str());
}
