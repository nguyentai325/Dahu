#ifndef DISPATCHER_HPP
# define DISPATCHER_HPP

# include <utility>
# include <vector>
# include "qattribute.hpp"
# include <mln/qt/imageviewer.hpp>
# include "plotwindow.hpp"

namespace mln
{

  class QDispatcher : public QObject
  {
    Q_OBJECT;

  public:
    template <typename V>
    QDispatcher(const image2d<V>& K,
		const image2d<unsigned>& parent,
		const std::vector<unsigned>& S);

    void addImageWindow(qt::ImageViewer* win);
    void addImageWindowToFilter(qt::ImageViewer* win, const image2d<rgb8>& mean);
    //void addAttribute(QAttributeBase* attr);
    void setPlotWindow(PlotWindow* pltwin);

  protected slots:
    void onNodeSelected(const mln::point2d& p);
    void onNodeSelected(const mln::image2d<bool>& mask);
    void onPointSelected(const mln::point2d& p);

  private:
    void doNodeSection(qt::ImageViewer*);
    void doFiltering(std::pair<qt::ImageViewer*, image2d<rgb8> >& obj);

    const image2d<unsigned>&	 m_parent;
    const std::vector<unsigned>& m_S;
    std::vector<point2d>	 m_leaves;
    bool			 m_leaf_attach;

    image2d<bool>		 m_mask_selection;
    //std::vector<QAttributeBase*> m_attributes;
    PlotWindow*			 m_pltwin;
    std::vector<qt::ImageViewer*> m_windows;
    std::vector< std::pair<qt::ImageViewer*, image2d<rgb8> > > m_fwins;
  };


  template <typename V>
  QDispatcher::QDispatcher(const image2d<V>& K,
			   const image2d<unsigned>& parent,
			   const std::vector<unsigned>& S)
    : QObject(),
      m_parent(parent),
      m_S (S)
  {
    resize(m_mask_selection, parent).init(true);
    m_leaves.reserve(m_S.size());

    {
      for (int i = S.size() - 1; i > 0; --i)
	{
	  unsigned x = S[i];
	  if (K[x] != K[parent[x]])
	    {
	      if (m_mask_selection[x])
		m_leaves.push_back(parent.point_at_index(x));
	      else
		m_mask_selection[parent[x]] = false;
	    }
	}
    }

  }

}

#endif // ! DISPATCHER_HPP
