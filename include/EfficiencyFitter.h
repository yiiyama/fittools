//-----------------------------------------------------------
// Unbinned maximum likelihood fit to pass / fail data
//-----------------------------------------------------------

#ifndef EfficiencyFitter_H
#define EfficiencyFitter_H

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#include <vector>
#include <utility>

class EfficiencyFitter {
 public:
  ~EfficiencyFitter();

  void setTarget(TTree*, char const* passdef, char const* xexpr, char const* yexpr = 0);
  void setLikelihood(RooAbsReal const* likelihood, RooArgList const* params, RooRealVar* xvar, RooRealVar* yvar = 0);

  int fit(int verbosity = 0);

  RooRealVar* getParam(char const*) const;
  double getNLL() const;
  void reset();

  void setXRange(double, double);
  void setYRange(double, double);

  bool hesse{false};
  bool minos{false};

  static EfficiencyFitter* singleton();
  static void FCN(int&, double*, double&, double*, int);

 private:
  EfficiencyFitter() {}
  double fcn(double const*) const;

  std::vector<RooRealVar*> parameters_{};
  RooAbsReal const* likelihood_{0};

  RooRealVar* xvar_{0};
  RooRealVar* yvar_{0};

  double* xmin_{0};
  double* xmax_{0};
  double* ymin_{0};
  double* ymax_{0};

  std::vector<std::pair<double, double>> pass_{};
  std::vector<std::pair<double, double>> fail_{};
};

#endif
