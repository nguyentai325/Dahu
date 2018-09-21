#ifndef APPS_G2_COMPUTE_CTOS_HPP
# define APPS_G2_COMPUTE_CTOS_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/core/colors.hpp>
# include <mln/core/grays.hpp>
# include <mln/morpho/component_tree/component_tree.hpp>

namespace mln
{

  enum e_ctos_attribute {
    CTOS_DEPTH,         // Depth of longest path to A
    CTOS_COUNT,         // Number of nodes including A
    CTOS_PW_COUNT,      // Number of nodes including x
  };

  struct ctos_extra_params_t
  {
    bool        export_marginal_depth = false;
    std::string export_marginal_depth_path;
  };


  /// \brief A convenient wrapper around:
  /// * The Graph of shapes computation + depth attribute valuation
  /// * The Saturated Maxtree algorithm
  ///
  /// \param input The input image
  /// \param[out] depth if not NULL, store the depth image inside.
  template <class V>
  morpho::component_tree<unsigned, image2d<unsigned> >
  compute_ctos(const image2d<V>& input,
               image2d<uint16>* depth = nullptr,
               e_ctos_attribute attribute = CTOS_DEPTH,
               ctos_extra_params_t params = ctos_extra_params_t ()
               );


  morpho::component_tree<unsigned, image2d<unsigned> >
  compute_ctos_from_maxtrees(const image2d<rgb8>& input,
                             image2d<uint16>* imdepth,
                             bool mintree = false);

  extern template
  morpho::component_tree<unsigned, image2d<unsigned> >
  compute_ctos<rgb8>(const image2d<rgb8>& input,
                     image2d<uint16>* depth,
                     e_ctos_attribute attribute,
                     ctos_extra_params_t params
                     );

  extern template
  morpho::component_tree<unsigned, image2d<unsigned> >
  compute_ctos<rgb16>(const image2d<rgb16>& input,
                      image2d<uint16>* depth,
                      e_ctos_attribute attribute,
                      ctos_extra_params_t params
                      );

}

#endif // ! APPS_G2_COMPUTE_CTOS_HPP
