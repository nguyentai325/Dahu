#include <QGridLayout>
#include <QMenuBar>
#include "plotwindow.hpp"

PlotWindow::PlotWindow()
{
  QWidget* main = new QWidget();

  this->init();

  main->setLayout(m_layout);
  this->setCentralWidget(main);
}

void
PlotWindow::register_attribute(Attribute* attribute)
{
  for (QString s : attribute->names())
    {
      m_attributes.push_back(attribute);
      m_attribute_selector->addItem(s);
    }
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
PlotWindow::init()
{
  m_layout = new QVBoxLayout(this);

  QVBoxLayout* m_panel_left = new QVBoxLayout;
  m_attribute_selector = new QComboBox(this);
  m_attribute_runnew = new QCheckBox("Display in new plot", this);
  m_panel_left->addWidget(m_attribute_selector);
  m_panel_left->addWidget(m_attribute_runnew);

  QGridLayout* m_attribute_panel_layout = new QGridLayout;
  m_options_panel = new QFormLayout(this);
  m_options_btncpt = new QPushButton("Run", this);

  m_attribute_panel = new QGroupBox("Attribute properties");
  m_attribute_panel_layout->addLayout(m_panel_left, 0, 0);
  m_attribute_panel_layout->addLayout(m_options_panel, 0, 1);
  m_attribute_panel->setLayout(m_attribute_panel_layout);


  QObject::connect(m_attribute_selector, SIGNAL(currentIndexChanged(int)),
		   this, SLOT(onAttributeSelected(int)));

  m_options_btncpt->installEventFilter(this);

  m_layout->addWidget(m_attribute_panel);
  m_has_selected_point = false;

  this->createActions();
  this->createMenus();
}

void
PlotWindow::createMenus()
{
  m_display_menu = this->menuBar()->addMenu("Display");
  m_display_menu->addAction(m_action_rm_plot);
}

void
PlotWindow::createActions()
{
  m_action_rm_plot = new QAction("Remove last plot", this);
  QObject::connect(m_action_rm_plot, SIGNAL(triggered()),
		   this, SLOT(rm_last_plot()));
}

void
PlotWindow::displayPlot(int i, bool new_plot, bool force_run)
{
  Attribute* attr = m_attributes[i];
  displayOptions(attr->parameters());
  mln::QAttributeBase* plot = attr->getPlot(m_attribute_selector->itemText(i));

  if (plot == NULL or force_run) {
    attr->run();
    plot = attr->getPlot(m_attribute_selector->itemText(i));

    if (m_has_selected_point)
      plot->plotNode(m_current_selected_point);
  }

  // Display plot: 3 case
  // 1. if the plot is already displayed => do nothing
  // 2. else if not displayed and not new => replace the last one
  // 3. create a new one
  for (int i = 0; i < m_layout->count(); ++i)
    if (m_layout->itemAt(i)->widget() == plot)
      return;


  if (!new_plot and m_layout->count() > 1) { // remove last  plot
    QLayoutItem* old = m_layout->takeAt(m_layout->count() - 1);
    old->widget()->setVisible(false);
  }
  m_layout->addWidget(plot);
  plot->setVisible(true);
}



void
PlotWindow::onAttributeSelected(int i)
{
  this->displayPlot(i, m_attribute_runnew->isChecked() , false);
}


void
PlotWindow::plotNode(const mln::point2d& p)
{
  m_has_selected_point = true;
  m_current_selected_point = p;
  for (unsigned i = 0; i < m_attributes.size(); ++i)
    {
      QString s = m_attribute_selector->itemText(i);
      mln::QAttributeBase* plot = m_attributes[i]->getPlot(s);
      if (plot)
	plot->plotNode(p);
    }
}

void
PlotWindow::onNodeSelected(const mln::point2d& p)
{
  Q_EMIT(nodeSelected(p));
}

void
PlotWindow::onNodeSelected(const mln::image2d<bool>& pts)
{
  Q_EMIT(nodeSelected(pts));
}


void
PlotWindow::displayOptions(const QMap<QString, Attribute::Parameter>& params)
{
  while (!m_options_panel->isEmpty()) {
    QLayoutItem* item = m_options_panel->takeAt(0); // clear the panel
    if (item->widget())
      item->widget()->setVisible(false);
  }
  for (const auto& param: params) {
    std::cout << param.label->text().toStdString() << std::endl;
    if (param.obj->widget())
      m_options_panel->addRow(param.label, param.obj->widget());
    else
      m_options_panel->addRow(param.label, param.obj->layout());

    param.label->setVisible(true);
    QWidget* w;
    if ((w = param.obj->widget()) != NULL)
      w->setVisible(true);
  }
  if (!params.isEmpty())
    {
      m_options_panel->addRow(m_options_btncpt);
      m_options_btncpt->setVisible(true);
    }
  m_options_panel->invalidate();
  m_options_panel->update();
  m_attribute_panel->update();
}

void
PlotWindow::rm_last_plot()
{
  std::cout << "Fuck" << std::endl;
  if (m_layout->count() > 1) { // remove old plot
    QLayoutItem* old = m_layout->takeAt(m_layout->count() - 1);
    old->widget()->setVisible(false);
  }
}

bool
PlotWindow::eventFilter(QObject *obj, QEvent *event)
{
  // When recompute button is clicked
  // Recompute the attribute and update the curve
  if (obj == m_options_btncpt and event->type() == QEvent::MouseButtonRelease)
    {
      int i = m_attribute_selector->currentIndex();
      displayPlot(i, false, true);
      return false;
    }

  return QObject::eventFilter(obj, event);
}
