#ifndef PLOTWINDOW_HPP
# define PLOTWINDOW_HPP

# include <QMainWindow>
# include <QGroupBox>
# include <QVBoxLayout>
# include <QComboBox>
# include <QFormLayout>
# include <QPushButton>
# include <QAction>
# include <QMenu>
# include <QCheckBox>
# include <vector>
# include "attribute.hpp"

/// \brief define a the window class that holds the attribute plots
/// and allows to select attributes to plot.
class PlotWindow : public QMainWindow
{
  Q_OBJECT;
public:
  PlotWindow();

  void register_attribute(Attribute* attribute);

public slots:
  /// \brief Slot to actualize the plotted branch when
  /// a point is clicked
  void plotNode(const mln::point2d& pt);

signals:
  /// \brief Signal emited when a single node is selected
  void nodeSelected(const mln::point2d& pt);

  /// \brief Signal emited when several nodes are selected.
  /// The mask is a boolean image
  void nodeSelected(const mln::image2d<bool>& pts);

private slots:
  /// \brief Slot when a child has a node selected
  void onNodeSelected(const mln::point2d& pt);

  /// \brief Slot when a child has nodes selected
  /// The mask is a boolean image
  void onNodeSelected(const mln::image2d<bool>& pts);


private:
  /// \brief initialize the window, menu...
  void init();
  void createMenus();
  void createActions();
  void displayPlot(int idxAttribute,
		   bool new_plot = false,
		   bool force_run = false);

private slots:
  void onAttributeSelected(int i);
  void displayOptions(const QMap<QString, Attribute::Parameter>& params);
  void rm_last_plot();

protected:
  bool eventFilter(QObject *obj, QEvent *event);

private:
  QVBoxLayout* m_layout;
  QGroupBox*   m_attribute_panel;
  QComboBox*   m_attribute_selector;
  QFormLayout* m_options_panel;
  QPushButton* m_options_btncpt;
  QCheckBox*   m_attribute_runnew;

  QMenu*       m_display_menu;
  QAction*     m_action_rm_plot;

  mln::point2d      m_current_selected_point;
  bool	            m_has_selected_point;
private:
  std::vector<Attribute*>		m_attributes;
};

#endif // ! PLOTWINDOW_HPP
