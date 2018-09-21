#include "attribute.hpp"

void
Attribute::setSignals(mln::QAttributeBase* attribute)
{
  QObject::connect(attribute,
		   SIGNAL(nodeSelected(const mln::point2d&)),
		   this,
		   SLOT(onNodeSelected(const mln::point2d&)));

   QObject::connect(attribute,
		    SIGNAL(nodeSelected(const mln::image2d<bool>&)),
		    this,
		    SLOT(onNodeSelected(const mln::image2d<bool>&)));
}


void
  Attribute::onNodeSelected(const mln::point2d& p)
{
  Q_EMIT(nodeSelected(p));
}

void
Attribute::onNodeSelected(const mln::image2d<bool>& pts)
{
  Q_EMIT(nodeSelected(pts));
}

