#ifndef UTILS_HPP
# define UTILS_HPP

# include <apps/tos/topology.hpp>

namespace mln
{

  /// \brief Transform a ToS in order to keep only the 2 faces in the tree
  template <typename T>
  std::tuple< image2d<T>, image2d<unsigned>, std::vector<unsigned> >
  remove_01face_from_tos(const image2d<T>& K, const image2d<unsigned>& parent, const std::vector<unsigned>& s);



  /********************************/
  /*** Implementation           ***/
  /********************************/


  template <typename T>
  std::tuple< image2d<T>, image2d<unsigned>, std::vector<unsigned> >
  remove_01face_from_tos(const image2d<T>& K_, const image2d<unsigned>& parent_,
			 const std::vector<unsigned>& S_)
  {
    box2d dom = K_.domain();
    dom.pmax = (dom.pmax + 1) / 2;

    image2d<T> K(dom);
    image2d<unsigned> parent(dom);
    std::vector<unsigned> S;
    S.reserve(dom.size());


    copy(K_ | rng::filter(K_.domain(), K1::is_face_2), K);

    for (unsigned x: S_) {
      point2d p = K_.point_at_index(x);
      if (K1::is_face_2(p))
	S.push_back(K.index_of_point(p/2));
    }

    mln_foreach(const point2d& p, K.domain()) {
      point2d q = K_.point_at_index(parent_(2*p));
      mln_assertion(K.domain().has(q/2));
      parent(p) = K.index_of_point(q/2);
    }

    return std::make_tuple(std::move(K), std::move(parent), std::move(S));
  }


}

#endif // ! UTILS_HPP
