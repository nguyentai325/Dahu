#ifndef APPS_TOSGUI_ATTRIBUTES_MEANINGFULLNESS_HPP
# define APPS_TOSGUI_ATTRIBUTES_MEANINGFULLNESS_HPP


# include <mln/core/image/image2d.hpp>
# include <mln/morpho/filtering.hpp>
# include <apps/tosgui/attribute.hpp>
# include <apps/tosgui/qattribute.hpp>
# include <apps/attributes/meaningfullness.hpp>
# include <apps/saliency/extinction.hpp>

template <typename T, typename V>
class MeaningFullnessAttribute : public Attribute
{

public:

  MeaningFullnessAttribute(const mln::image2d<T>& ima_,
			   const mln::image2d<V>& K_,
			   const mln::image2d<unsigned>& parent_,
			   const std::vector<unsigned>& S_)
    : ima (ima_), K(K_), parent(parent_), S(S_),
      m_attribute_1 (NULL),
      m_attribute_2 (NULL),
      m_attribute_3 (NULL),
      m_attribute_4 (NULL),
      m_attribute_5 (NULL)
  {
    Parameter p1(new QLabel("epsilon"), new QLineEdit("5"));
    Parameter p2(new QLabel("alpha (Eint)"), new QLineEdit("0.3"));
    Parameter p3(new QLabel("beta (Econ)"), new QLineEdit("1"));
    Parameter p4(new QLabel("gamma (Econ)"), new QLineEdit("0.005"));

    m_params.insert(QString("eps"), p1);
    m_params.insert(QString("alpha"), p2);
    m_params.insert(QString("beta"), p3);
    m_params.insert(QString("gamma"), p4);
  }

  virtual
  ~MeaningFullnessAttribute()
  {
    if (m_attribute_1)
      {
	delete m_attribute_1;
	delete m_attribute_2;
	delete m_attribute_3;
	delete m_attribute_4;
	delete m_attribute_5;
      }
  }

  virtual
  std::vector<QString>
  names() const
  {
    std::vector<QString> nm;
    nm.push_back("MeaningFullness (Internal Energy Term - Variance)");
    nm.push_back("MeaningFullness (External Energy term - Curvature)");
    nm.push_back("MeaningFullness (Constraint Term - exp(-gamma.area) )");
    nm.push_back("MeaningFullness (All terms)");
    nm.push_back("MeaningFullness (Extinction)");
    return nm;
  }

  virtual
  QMap<QString, Parameter>&
  parameters()
  {
    return m_params;
  }

protected slots:

  virtual
  void showinfo(const mln::point2d& node)
  {
    float Vin = m_acc(node).cvar(m_acc(node).m_v_sum_int,
				 m_acc(node).m_v_sum_int_sqr,
				 m_acc(node).m_v_n_int);
    float Vout = m_acc(node).cvar(m_acc(node).m_v_sum_ext,
				  m_acc(node).m_v_sum_ext_sqr,
				  m_acc(node).m_v_n_ext);
    float Vtotal = m_acc(node).cvar(m_acc(node).m_v_sum_ext + m_acc(node).m_v_sum_int,
				  m_acc(node).m_v_sum_ext_sqr + m_acc(node).m_v_sum_int_sqr,
				  m_acc(node).m_v_n_ext + m_acc(node).m_v_n_int);

    auto c1 =  m_acc(node).m_v_sum_int / (float) m_acc(node).m_v_n_int;
    auto c2 =  m_acc(node).m_v_sum_ext / (float) m_acc(node).m_v_n_ext;
    auto c3 =  (m_acc(node).m_v_sum_ext + m_acc(node).m_v_sum_int) /
      (float) (m_acc(node).m_v_n_ext + m_acc(node).m_v_n_int);

    float a, b, delta;
    std::tie(a,b) = m_acc(node).inter(c1, Vin, c2, Vout, delta);

    std::cout << "Internal variance: " << Vin << std::endl
	      << "External variance: " << Vout << std::endl
	      << "Total variance: " << Vtotal << std::endl
	      << "Internal center : " << c1 << std::endl
	      << "External center : " << c2 << std::endl
	      << "Total center : " << c3 << std::endl
      ;
    std::cout << "a:" << a << "  b: " << b << std::endl;

    auto hf = [c1, Vin] (double x) { return 0.5 * (1 + std::erf( (x-c1)/std::sqrt(2*M_PI*Vin))); };
    auto hg = [c2, Vout] (double x) { return 0.5 * (1 + std::erf( (x-c2)/std::sqrt(2*M_PI*Vout))); };

    double if1, if2, if3, ig1, ig2, ig3;
      if1 = hf(a);
      ig1 = hg(a);
      if2 = hf(b) - if1;
      ig2 = hg(b) - ig1;
      if3 = 1 - if1 - if2;
      ig3 = 1 - ig1 - ig2;
      float m1, M1, m2, M2, m3, M3;
      std::tie(m1, M1) = std::minmax(if1, ig1);
      std::tie(m2, M2) = std::minmax(if2, ig2);
      std::tie(m3, M3) = std::minmax(if3, ig3);
      std::cout << "(" << m1 << "," << m2 << "," << m3 << ")" << std::endl;
      std::cout << "(" << M1 << "," << M2 << "," << M3 << ")" << std::endl;
  }


public:
  virtual
  void
  run()
  {
    int eps =  ((const QLineEdit*) (m_params["eps"].obj->widget()))->text().toInt();
    float alpha = ((const QLineEdit*) (m_params["alpha"].obj->widget()))->text().toFloat();
    float beta = ((const QLineEdit*) (m_params["beta"].obj->widget()))->text().toFloat();
    float gamma = ((const QLineEdit*) (m_params["gamma"].obj->widget()))->text().toFloat();

    if (m_attribute_1 == NULL)
      {
	m_eps = -eps;
	mln::resize(m_acc, K);
	mln::resize(m_con, K);
	mln::resize(m_curv, K);
	mln::resize(m_ext, K);

	m_attribute_1 = new mln::QAttribute<float>(m_con, parent, QString("Meaningfullness - Constraint: exp(-gamma.area)"));
	m_attribute_2 = new mln::QAttribute<float>(m_curv, parent, QString("Meaningfullness - Curvature"));
	m_attribute_3 = new mln::QAttribute<float>(m_ext, parent, QString("Meaningfullness - (VarExt + VarInt) / VarTot"));
	m_attribute_4 = new mln::QAttribute<float>(m_energy, parent, QString("Meaningfullness - energy"));
	m_attribute_5 = new mln::QAttribute<float>(m_extinction, parent, QString("Meaningfullness - extinction"));


	QObject::connect(m_attribute_3, SIGNAL(nodeSelected(const mln::point2d&)),
			 this, SLOT(showinfo(const mln::point2d&)));

	this->setSignals(m_attribute_1);
	this->setSignals(m_attribute_2);
	this->setSignals(m_attribute_3);
	this->setSignals(m_attribute_4);
	this->setSignals(m_attribute_5);
      }

    if (eps != m_eps) // recompute: eps has changed
      {
	m_energy = mln::meaningfullness(ima, K, parent, S, m_acc, alpha, beta, gamma, eps);
	m_eps = eps;

	// Compute aux data
	{
	  mln::transform(m_acc,
			 [beta, gamma] (const mln::internal::energy_t<T>& aux) { return beta * std::exp(-gamma * aux.m_area); },
			 m_con);

	  mln::transform(m_acc,
			 [alpha] (const mln::internal::energy_t<T>& aux) { return (float) alpha * aux.m_e_sumcurv / aux.m_e_length; },
			 m_curv);

	  mln::transform(m_acc,
			 [] (const mln::internal::energy_t<T>& aux) { return aux.external_energy(); },
			 m_ext);
	}
	m_extinction = extinction(m_energy, K, parent, S, std::less<float>());
      }
    else
      {
	// We juste have new parameters for alpha and beta
	  mln::transform(m_acc,
			 [beta, gamma] (const mln::internal::energy_t<T>& aux) { return beta * std::exp(-gamma * aux.m_area); },
			 m_con);

	  mln::transform(m_acc,
			 [alpha] (const mln::internal::energy_t<T>& aux) { return (float) alpha * aux.m_e_sumcurv / aux.m_e_length; },
			 m_curv);

	  mln::transform(m_acc,
			 [] (const mln::internal::energy_t<T>& aux) { return aux.external_energy(); },
			 m_ext);

	  std::cout << "a: " << alpha << " a0: " << beta << " a1: " << gamma << std::endl;
	  mln::transform(m_acc,
			 [alpha, beta, gamma] (const mln::internal::energy_t<T>& aux) { return aux.energy(alpha, beta, gamma); },
			 m_energy);
	  m_extinction = extinction(m_energy, K, parent, S, std::less<float>());
      }
  }

  virtual
  mln::QAttributeBase*
  getPlot(const QString& name)
  {
    static const char* names[] = {
      "MeaningFullness (Internal Energy Term - Variance)",
      "MeaningFullness (External Energy term - Curvature)",
      "MeaningFullness (Constraint Term - exp(-gamma.area) )",
      "MeaningFullness (All terms)",
      "MeaningFullness (Extinction)"
    };

    if (name == names[0])
      return m_attribute_3;
    else if (name == names[1])
      return m_attribute_2;
    else if (name == names[2])
      return m_attribute_1;
    else if (name == names[3])
      return m_attribute_4;
    else
      return m_attribute_5;
  }


private:
  const mln::image2d<T>&        ima;
  const mln::image2d<V>&        K;
  const mln::image2d<unsigned>& parent;
  const std::vector<unsigned>&  S;

  QMap<QString, Parameter>	m_params;
  int	m_eps;

  mln::image2d< mln::internal::energy_t<T> > m_acc;
  mln::image2d<float>		m_con;
  mln::image2d<float>		m_curv;
  mln::image2d<float>		m_ext;
  mln::image2d<float>		m_energy;
  mln::image2d<float>		m_extinction;

  mln::QAttribute<float>*	m_attribute_1;
  mln::QAttribute<float>*	m_attribute_2;
  mln::QAttribute<float>*	m_attribute_3;
  mln::QAttribute<float>*	m_attribute_4;
  mln::QAttribute<float>*	m_attribute_5;
};

#endif // !APPS_TOSGUI_ATTRIBUTES_MEANINGFULLNESS_HPP
