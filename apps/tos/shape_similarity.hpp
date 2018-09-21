#ifndef SHAPE_SIMILARITY_HPP
# define SHAPE_SIMILARITY_HPP

# include <vector>
# include <mln/core/image/image2d.hpp>
# include <mln/morpho/saturate.hpp>

namespace mln
{

  template <typename L, typename V>
  float
  coveryRate(const image2d<L>& labels, int nlabels,
	     const image2d<V>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& S)
  {
    mln_precondition(labels.domain() == K.domain());

    image2d< std::vector<int> > inter;
    resize(inter, K);

    inter[S[0]].resize(nlabels + 1);
    for (unsigned x: S)
      if (K[x] != K[parent[x]])
	inter[x].resize(nlabels + 1);

    // Number of pixel per zone
    std::vector<unsigned> count(nlabels+1);

    // Count for each shape the number of pixels of each label i.e. Card(Shape Inter Zone)
    for (int l = 0; l <= nlabels; ++l)
      {
	auto sat = morpho::saturate(labels == l, c4, point2d{0,0});

	mln_foreach(auto& px, sat.pixels())
	  if (px.val())
	    {
	      ++count[l];
	      unsigned k = K.index_of_point(px.point());
	      unsigned x = K[k] == K[parent[k]] ? parent[k] : k;
	      while (x != parent[x]) {
		++inter[x][l];
		x = parent[x];
	      }
	      ++inter[x][l]; // root
	    }
      }

    // Compute the cardinality of each shape
    image2d<unsigned> area;
    {
      resize(area, K);
      area[S[0]] = 1;
      for (int i = S.size()-1; i > 0; --i)
	{
	  unsigned x = S[i];
	  area[x] += 1;
	  area[parent[x]] += area[x];
	}
    }


    // Compute the Card(A Inter B) / Card(A U B)
    std::vector<float> covery(nlabels+1, 0);
    {
      for (unsigned x: S)
	{
	  if (K[x] != K[parent[x]] or x == parent[x]) {
	    //std::cout << "Shape: " << x << ";" << K[x] << " size: " << area[x] << std::endl;
	    for (int i = 0; i <= nlabels; ++i) {
	      float cov = (float) inter[x][i] / (area[x] + count[i] - inter[x][i]);
	      covery[i] = std::max(covery[i], cov);
	      //std::cout << "\t Region: " << i << " size: " << count[i] << " inter: " << inter[x][i] << " cov: " << cov << std::endl;
	    }
	  }
	}
    }

    // display: covery
    {
      for (int i = 0; i <= nlabels; ++i)
	std::cout << "Region " << i << " : " << covery[i] << std::endl;
    }

    return std::accumulate(covery.begin(), covery.end(), 0.0);
  }

}

#endif // ! SHAPE_SIMILARITY_HPP
