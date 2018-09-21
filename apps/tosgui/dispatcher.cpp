#include <algorithm>

#include "dispatcher.hpp"
#include <mln/core/algorithm/fill.hpp>
#include <mln/core/image/sub_image.hpp>

namespace mln
{

  void
  QDispatcher::addImageWindow(qt::ImageViewer* win)
  {
    m_windows.push_back(win);
    QObject::connect(win, SIGNAL(pointSelected(const mln::point2d&)),
		     this, SLOT(onPointSelected(const mln::point2d&)));
  }

  void
    QDispatcher::addImageWindowToFilter(qt::ImageViewer* win,
                                        const image2d<rgb8>& mean)
  {
    m_fwins.push_back(std::make_pair(win, mean));
    QObject::connect(win, SIGNAL(pointSelected(const mln::point2d&)),
		     this, SLOT(onPointSelected(const mln::point2d&)));
  }

  // void
  // QDispatcher::addAttribute(QAttributeBase* attr)
  // {
  //   m_attributes.push_back(attr);

  //   QObject::connect(attr, SIGNAL(nodeSelected(const point2d&)),
  // 		     this, SLOT(onNodeSelected(const point2d&)));
  //   QObject::connect(attr, SIGNAL(nodeSelected(const image2d<bool>&)),
  // 		     this, SLOT(onNodeSelected(const image2d<bool>&)));
  // }

  void
  QDispatcher::setPlotWindow(PlotWindow* pltwin)
  {
    m_pltwin = pltwin;
    QObject::connect(pltwin, SIGNAL(nodeSelected(const mln::point2d&)),
		     this, SLOT(onNodeSelected(const mln::point2d&)));
    QObject::connect(pltwin, SIGNAL(nodeSelected(const mln::image2d<bool>&)),
		     this, SLOT(onNodeSelected(const mln::image2d<bool>&)));
  }



  void
  QDispatcher::onNodeSelected(const point2d& p)
  {
    unsigned x = m_parent.index_of_point(p);

    fill(m_mask_selection, false);
    for (unsigned p: m_S)
      if (p == x or m_mask_selection[m_parent[p]])
	m_mask_selection[p] = true;

    for (qt::ImageViewer* win : m_windows)
      this->doNodeSection(win);

    for (auto& win : m_fwins)
      this->doFiltering(win);
  }

  void
  QDispatcher::onNodeSelected(const image2d<bool>& mask)
  {
    fill(m_mask_selection, false);
    int n = 0;
    for (unsigned p: m_S)
      if (mask[p] or m_mask_selection[m_parent[p]])
	{
	  m_mask_selection[p] = true;
	  ++n;
	}

    std::cout << "Filtering: " << n << std::endl;

    for (qt::ImageViewer* win : m_windows)
      this->doNodeSection(win);

    for (auto& win : m_fwins)
      this->doFiltering(win);
  }


  void
  QDispatcher::onPointSelected(const point2d& p)
  {
    // Select leaf
    auto dist = [&p] (const point2d& q) {
      return std::abs(p[0] - q[0]) + std::abs(p[1] - q[1]);
    };

    std::vector<unsigned> tmp(m_leaves.size());
    std::transform(m_leaves.begin(), m_leaves.end(), tmp.begin(), dist);
    auto x = std::min_element(tmp.begin(), tmp.end());

    point2d pmin = m_leaves[x - tmp.begin()];

    // Transmit selection
    if (m_pltwin)
      m_pltwin->plotNode(pmin);
  }

  void
  QDispatcher::doNodeSection(qt::ImageViewer* win)
  {
    image2d<rgb8>& view = win->getView();
    win->reset();

    auto x = view | where(m_mask_selection);
    mln_foreach(auto& v, x.values()) {
      v[0] = v[0] / 2 + 128;
    }

    win->update();
  }

  void
  QDispatcher::doFiltering(std::pair<qt::ImageViewer*, image2d<rgb8> >& obj)
  {
    qt::ImageViewer* win = obj.first;
    image2d<rgb8>& view = obj.first->getView();
    image2d<rgb8>& mean = obj.second;
    win->reset();

    for (unsigned x: m_S)
      if (m_mask_selection[m_parent[x]]) {
        view[x] = view[m_parent[x]];
      } else if (m_mask_selection[x]) {
        view[x] = mean[m_parent[x]];
      }

    win->update();
  }

}
