#ifndef APPS_TOSGUI_ATTRIBUTE_HPP
# define APPS_TOSGUI_ATTRIBUTE_HPP

# include <QString>
# include <QLabel>
# include <QWidget>
# include <QLayoutItem>
# include <QMap>
# include <vector>
# include "qattribute.hpp"

/// \brief Abstract class for attributes
class Attribute : public QObject
{
  Q_OBJECT;
public:
  struct Parameter
  {
    enum etype { VALUE = 0, CHOICE = 1};

    Parameter() = default;

    Parameter(QLabel* lbl, QWidget* w)
      : label (lbl), obj (new QWidgetItem(w))
    {
    }

    Parameter(QLabel* lbl, QLayoutItem* w)
      : label (lbl), obj (w)
    {
    }

    QLabel*  label;
    QLayoutItem* obj;
    //etype    type;

    // ~Parameter()
    // {
    //   delete label;
    //   delete obj;
    // }
  };

  /// \brief return the attribute's names
  virtual std::vector<QString> names() const = 0;

  /// \brief return the parameters of the attributes (if so)
  virtual QMap<QString, Parameter>& parameters() = 0;

  /// \brief runs the method with teh following parameters
  virtual void run() = 0;

  /// \brief returns the curve of the attribute
  virtual mln::QAttributeBase* getPlot(const QString& attribute_name) = 0;

signals:
  void nodeSelected(const mln::point2d&);
  void nodeSelected(const mln::image2d<bool>&);


protected:
  void setSignals(mln::QAttributeBase* attribute);

protected slots:
  void onNodeSelected(const mln::point2d&);
  void onNodeSelected(const mln::image2d<bool>&);

  virtual
  void showinfo(const mln::point2d&)
  {
  }

};

#endif // ! APPS_TOSGUI_ATTRIBUTE_HPP
