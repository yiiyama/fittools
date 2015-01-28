#ifndef IntegralChi2Fitter_H
#define IntegralChi2Fitter_H

#include "TH1.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include <vector>
#include <utility>

class IntegralChi2Fitter {
 public:
  ~IntegralChi2Fitter();

  void setTarget(TH1 const*, bool = true, char const* = 0);
  void setTarget(int, double const*, double const*, double const*, bool = true, char const* = 0);
  void setTarget(int, double, double, double const*, double const*, bool = true, char const* = 0);
  void setFunction(RooAbsReal const*, RooRealVar*, char const* = 0, RooArgSet* = 0);
  void setFunction(RooAbsReal const*, RooArgList const*, char const* = 0, RooArgSet* = 0);
  int fit(int = 0);

  RooRealVar* getParam(char const*) const;
  double getChi2() const;
  TMatrixDSym const* getCovarianceMatrix() const { return covariance_; }
  unsigned getNdof() const;
  void excludeBin(int);
  void reset();

  TH1* histogram() const;
  TH1* pullHistogram() const;
  double integral() const;
  double sumTarget() const;

  bool hesse;
  bool minos;

  static IntegralChi2Fitter* singleton();
  static void FCN(int&, double*, double&, double*, int);

 private:
  struct Bin {
    Bin();
    Bin(unsigned);
    Bin(Bin const&);
    ~Bin();
    Bin& operator=(Bin const&);

    void setIndex(unsigned, int);
    int index(unsigned) const;

    unsigned ndim;
    TString name;
    double content;
    double error;
    double volume;
    RooAbsReal* integral;

  private:
    unsigned* index_;
  };

  IntegralChi2Fitter();
  void clearFunction();
  double fcn(double const*) const;

  int dimension_;
  TString functionName_;

  std::vector<int> binIndices_;

  std::vector<Bin> bins_;

  std::vector<std::vector<double> > binEdges_;

  std::vector<RooRealVar*> params_;

  TMatrixDSym* covariance_;

  RooArgSet observables_;
  RooAbsPdf const* extendPdf_;
};

#endif
