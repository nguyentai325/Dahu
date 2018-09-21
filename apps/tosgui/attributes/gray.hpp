#ifndef APPS_TOSGUI_ATTRIBUTES_GRAY_HPP
# define APPS_TOSGUI_ATTRIBUTES_GRAY_HPP

# include <mln/core/image/image2d.hpp>
# include <apps/tosgui/attribute.hpp>
# include <apps/tosgui/qattribute.hpp>

template <typename V>
class GrayLevelAttribute : public Attribute
{
  static_assert(std::is_convertible<V, float>::value, "V must convertible to float");

public:

  GrayLevelAttribute(const mln::image2d<V>& K_,
		     const mln::image2d<unsigned>& parent_,
		     const std::vector<unsigned>& S_)
    : K(K_), parent(parent_), S(S_),
      m_attribute (NULL)
  {
  }

  virtual
  ~GrayLevelAttribute()
  {
    if (m_attribute)
      delete m_attribute;
  }

  virtual
  std::vector<QString>
  names() const
  {
    std::vector<QString> nm;
    nm.push_back("Gray Level");
    return nm;
  }

  virtual
  QMap<QString, Parameter>&
  parameters()
  {
    static QMap<QString, Parameter> params;
    return params;
  }

  virtual
  void
  run()
  {
    if (m_attribute == NULL) {
      m_attribute = new mln::QAttribute<V>(K, parent, QString("Gray Levels"));
      this->setSignals(m_attribute);
    }
  }

  virtual
  mln::QAttributeBase*
  getPlot(const QString&)
  {
    return m_attribute;
  }


private:
  const mln::image2d<V>&        K;
  const mln::image2d<unsigned>& parent;
  const std::vector<unsigned>&  S;

  mln::QAttribute<V>*		m_attribute;
};

#endif // ! APPS_TOSGUI_ATTRIBUTES_GRAY_HPP
