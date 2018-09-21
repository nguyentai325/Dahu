#ifndef MLN_QT_QTIMAGE_HPP
# define MLN_QT_QTIMAGE_HPP

# include <QImage>
# include <QColor>
# include <mln/core/image/image2d.hpp>
# include <mln/core/grays.hpp>
# include <mln/core/colors.hpp>
# include <mln/core/algorithm/transform.hpp>
# include <mln/core/algorithm/copy.hpp>

namespace mln
{
  namespace qt
  {

    class QtImageBase
    {
    public:
      QtImageBase(int nrows, int ncols, int border);
      virtual ~QtImageBase() = default;

      const image2d<rgb8>&	getView() const;
      image2d<rgb8>&		getView();
      const QImage&		getQImage() const;
      QImage&			getQImage();

      virtual void reset() = 0;


      void update();

    protected:
      image2d<rgb8>	m_view;
      QImage		m_qima;

      static QVector<QRgb> default_lut8;
    };


    template <typename V>
    class QtImage : public QtImageBase
    {
    public:
      QtImage(const image2d<V>& ima);


      virtual void reset();

    private:
      const image2d<V>& m_ima;
    };


    /******************************************/
    /****          Implementation          ****/
    /******************************************/

    namespace internal
    {
      // FIXME overload or specialization ?

      template <typename V>
      struct qt_format
      {
        static const QImage::Format format = QImage::Format_Invalid;
      };

      template <>
      struct qt_format<uint8>
      {
        static const QImage::Format format = QImage::Format_Indexed8;
      };

      template <>
      struct qt_format<rgb8>
      {
        static const QImage::Format format = QImage::Format_RGB888;
      };

      template <typename V>
      void rgb_convert(const image2d<V>& ima, image2d<rgb8>& out);

      template <>
      inline
      void rgb_convert<uint8>(const image2d<uint8>& ima, image2d<rgb8>& out)
      {
	transform(ima, [] (const uint8& x) { return rgb8{x,x,x}; }, out);
      }

      template <>
      inline
      void rgb_convert<rgb8>(const image2d<rgb8>& ima, image2d<rgb8>& out)
      {
	copy(ima, out);
      }

    }


    template <typename V>
    QtImage<V>::QtImage(const image2d<V>& ima)
      : QtImageBase(ima.nrows(), ima.ncols(), ima.border()),
        m_ima(ima)
    {
      reindex(this->m_view, m_ima);

      if (internal::qt_format<V>::format == QImage::Format_Indexed8)
        {
          if (this->default_lut8.empty()) {
            this->default_lut8.resize(256);
            for (int i = 0; i < 256; ++i)
              this->default_lut8[i] = QColor(i,i,i).rgb();
          }
          //this->setColorTable(this->default_lut8);
        }
      this->reset();
    }

    template <typename V>
    void
    QtImage<V>::reset()
    {
      internal::rgb_convert(m_ima, this->m_view);
    }

  } // end of namespace mln::qt

} // end of namespace mln

#endif //!MLN_QT_QTIMAGE_HPP
