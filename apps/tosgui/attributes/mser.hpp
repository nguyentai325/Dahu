#ifndef APPS_TOSGUI_ATTRIBUTES_MSER_HPP
# define APPS_TOSGUI_ATTRIBUTES_MSER_HPP

# include <mln/core/image/image2d.hpp>
# include <mln/morpho/filtering.hpp>
# include <apps/tos/topology.hpp>
# include <apps/tosgui/attribute.hpp>
# include <apps/tosgui/qattribute.hpp>
# include <apps/attributes/MSER.hpp>
# include <apps/saliency/closure.hpp>

# include <QComboBox>
# include <QVBoxLayout>
# include <QLabel>
# include <QComboBox>

template <typename V>
class MSERAttribute : public Attribute
{

public:

  MSERAttribute(const mln::image2d<V>& K_,
		     const mln::image2d<unsigned>& parent_,
		     const std::vector<unsigned>& S_)
    : K(K_), parent(parent_), S(S_),
      m_param_h (-1),
      m_attribute_1 (NULL),
      m_attribute_2 (NULL),
      m_attribute_3 (NULL),
      m_attribute_4 (NULL)
  {
    Parameter p1(new QLabel("h"), new QLineEdit("20"));

    QComboBox* p2cb = new QComboBox;
    p2cb->addItem("MSER_DIFF: A(par) - A(cur)");
    p2cb->addItem("MSER_RATIO: (A(CUR) / A(par)" );
    p2cb->addItem("MSER_NORM: (A(par) / A(cur) - 1)");
    p2cb->setCurrentIndex(1);

    Parameter p2(new QLabel("MSER mode"), p2cb);
    Parameter p3(new QLabel("a0"), new QLineEdit("5"));
    Parameter p4(new QLabel("a1"), new QLineEdit("0.0005"));

    m_params.insert(QString("h"), p1);
    m_params.insert(QString("MSER mode"), p2);
    m_params.insert(QString("a0"), p3);
    m_params.insert(QString("a1"), p4);

    // Closing form
    QVBoxLayout* b = new QVBoxLayout;

    QComboBox* ch = new QComboBox;
    ch->addItem("Area");
    ch->addItem("Height");

    b->addWidget(new QCheckBox("Active"));
    b->addWidget(ch);
    b->addWidget(new QLineEdit("10"));

    Parameter p5(new QLabel("Closing (Shape Space)"), b);
    m_params.insert(QString("closing"), p5);
  }

  virtual
  ~MSERAttribute()
  {
    if (m_attribute_1)
      {
	delete m_attribute_1;
	delete m_attribute_2;
	delete m_attribute_3;
	delete m_attribute_4;
      }
  }

  virtual
  std::vector<QString>
  names() const
  {
    std::vector<QString> nm;
    nm.push_back("MSER (Delta area term)");
    nm.push_back("MSER (Constraint term)");
    nm.push_back("MSER (Delta + Constraint term)");
    nm.push_back("MSER (Closure of Energy)");
    return nm;
  }

  virtual
  QMap<QString, Parameter>&
  parameters()
  {
    return m_params;
  }

  virtual
  void
  run()
  {
    if (m_area.domain().size() == 0)
      m_area = mln::morpho::area_compute(K, parent, S, mln::K1::is_face_2);

    if (m_attribute_1 == NULL) {
      m_attribute_1 = new mln::QAttribute<float>(m_mser, parent, QString("MSER (Delta area)"));
      m_attribute_2 = new mln::QAttribute<float>(m_cons, parent, QString("MSER (Constraint term)"));
      m_attribute_3 = new mln::QAttribute<float>(m_mser_, parent, QString("MSER (Delta + Constraint term)"));
      m_attribute_4 = new mln::QAttribute<float>(m_clo, parent, QString("MSER (Closure of Energy)"));

      this->setSignals(m_attribute_1);
      this->setSignals(m_attribute_2);
      this->setSignals(m_attribute_3);
      this->setSignals(m_attribute_4);
    }

    // const QLineEdit* htxt = (QLineEdit*) (params["h"].obj);
    // std::cout << htxt << std::endl;
    // std::cout << htxt->text().toStdString() << std::endl;
    int h =  ((const QLineEdit*) (m_params["h"].obj->widget()))->text().toInt();
    eMSER_attribute mode = (eMSER_attribute) ((const QComboBox*) (m_params["MSER mode"].obj->widget()))->currentIndex();
    float a0 = ((const QLineEdit*) (m_params["a0"].obj->widget()))->text().toFloat();
    float a1 = ((const QLineEdit*) (m_params["a1"].obj->widget()))->text().toFloat();

    std::cout << "Params: " << h << " / " << mode << std::endl;

    if (h != m_param_h)
      {
	m_mser     = compute_MSER_attribute(K, K, parent, S, h, mode);
	m_cons     = mln::transform(m_area, [a0,a1](float x) -> float {
	    return a0 * std::exp(-x * a1);
	  });
	m_mser_    = mln::eval(mln::transform(mln::imzip(m_mser, m_cons),
					      [](const std::tuple<float,float>& v) {
						return std::max(std::get<0>(v), std::get<1>(v));
					      }));
      }

    std::cout << "Fuck" << std::endl;
    QLayout* l = m_params["closing"].obj->layout();
    bool param_clo_active = ((const QCheckBox*) (l->itemAt(0)->widget()))->isChecked();
    int  param_clo_attribute = ((const QComboBox*) (l->itemAt(1)->widget()))->currentIndex();
    float param_clo_value = ((const QLineEdit*) (l->itemAt(2)->widget()))->text().toFloat();
    if (param_clo_active)
      {
	if (param_clo_attribute == 0)
	  m_clo = area_close(m_mser_, K, parent, S, param_clo_value);
	else
	  m_clo = height_close(m_mser_, K, parent, S, param_clo_value);
      }
    else
      m_clo = m_mser_;

    m_param_h = h;
  }

  virtual
  mln::QAttributeBase*
  getPlot(const QString& name)
  {
    if (name == "MSER (Delta area term)")
      return m_attribute_1;
    else if (name == "MSER (Constraint term)")
      return m_attribute_2;
    else if (name == "MSER (Delta + Constraint term)")
      return m_attribute_3;
    else
      return m_attribute_4;
  }


private:
  const mln::image2d<V>&        K;
  const mln::image2d<unsigned>& parent;
  const std::vector<unsigned>&  S;

  int				m_param_h;

  QMap<QString, Parameter>	m_params;
  mln::image2d<unsigned>	m_area;
  mln::image2d<float>		m_mser;
  mln::image2d<float>		m_cons;
  mln::image2d<float>		m_mser_;
  mln::image2d<float>		m_clo;

  mln::QAttribute<float>*	m_attribute_1;
  mln::QAttribute<float>*	m_attribute_2;
  mln::QAttribute<float>*	m_attribute_3;
  mln::QAttribute<float>*	m_attribute_4;
};


#endif
