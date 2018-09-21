#include <mln/qt/qtimage.hpp>


namespace mln
{

  namespace qt
  {

    QtImageBase::QtImageBase(int nrows, int ncols, int border)
      : m_view(nrows, ncols, border),
	m_qima((uchar*)&m_view(m_view.domain().pmin),
	       ncols,
	       nrows,
	       m_view.strides()[0],
	       QImage::Format_RGB888)
    {
    }

    const image2d<rgb8>&
    QtImageBase::getView() const
    {
      return m_view;
    }

    image2d<rgb8>&
    QtImageBase::getView()
    {
      return m_view;
    }

    const QImage&
    QtImageBase::getQImage() const
    {
      return m_qima;
    }

    QImage&
    QtImageBase::getQImage()
    {
      return m_qima;
    }

    void
    QtImageBase::update()
    {
      m_qima = QImage((uchar*)&m_view(m_view.domain().pmin),
                      m_view.ncols(),
                      m_view.nrows(),
                      m_view.strides()[0],
                      QImage::Format_RGB888);
    }

    QVector<QRgb> QtImageBase::default_lut8;
  }

}
