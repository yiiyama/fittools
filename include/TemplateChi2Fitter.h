#ifndef TemplateChi2Fitter_H
#define TemplateChi2Fitter_H

#include "TH1.h"
#include "TString.h"
#include "TObjArray.h"
#include "TDirectory.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

#include <vector>

class TemplateChi2Fitter {
 public:
  ~TemplateChi2Fitter();

  void setTarget(TH1 const*);
  void setTarget(TObjArray const*);
  void setTarget(int, double const*, double const*);
  void addTemplate(TH1 const*, char const* = 0, RooAbsReal* = 0);
  void addTemplate(TObjArray const*, char const* = 0, RooAbsReal* = 0);
  void addTemplate(int, double const*, double const*, char const* = 0, RooAbsReal* = 0);
  void changeTarget(TH1 const*, unsigned = 0);
  int fit(int = 0);

  RooRealVar* getFloat(char const*) const;
  double getChi2() const;
  unsigned getNdof() const;
  void plot(TDirectory*) const;
  void excludeBin(int, int = 1, int = 0);

  static TemplateChi2Fitter* singleton();
  static void FCN(int&, double*, double&, double*, int);

 private:
  TemplateChi2Fitter();
  double fcn(double const*) const;

  std::vector<int> nX_;
  std::vector<int> nY_;
  std::vector<double> targetContents_;
  std::vector<double> targetErrors_;
  std::vector<std::vector<double> > templateContents_;
  std::vector<std::vector<double> > templateErrors_;
  std::vector<unsigned> binOffsets_;
  std::vector<unsigned> bins_;
  RooRealVar chi2_;
  std::vector<RooRealProxy*> scales_;
  std::vector<RooRealVar*> floats_;
};

#endif
