#include "../include/EfficiencyFitter.h"

#include "TTree.h"
#include "TVirtualFitter.h"
#include "RooArgSet.h"

#include <algorithm>
#include <cmath>
#include <cstring>

EfficiencyFitter::~EfficiencyFitter()
{
  reset();
}

void
EfficiencyFitter::setTarget(TTree* _tree, char const* _passdef, char const* _xexpr, char const* _yexpr/* = 0*/)
{
  if (!_tree)
    return;

  bool hasY(_yexpr && std::strlen(_yexpr) != 0);
  
  TString expr(_passdef);
  expr += " != 0:";
  expr += _xexpr;
  if (hasY) {
    expr += ":";
    expr += _yexpr;
  }

  pass_.clear();
  fail_.clear();

  _tree->SetEstimate(_tree->GetEntries() + 1);

  long nEntries(_tree->Draw(expr, "", "goff"));

  double* passArr(_tree->GetV1());
  double* xArr(_tree->GetV2());
  double* yArr(_tree->GetV3());

  for (long iE(0); iE != nEntries; ++iE) {
    if (xmin_ && xArr[iE] < *xmin_)
      continue;
    if (xmax_ && xArr[iE] > *xmax_)
      continue;
    if (hasY) {
      if (ymin_ && yArr[iE] < *ymin_)
        continue;
      if (ymax_ && yArr[iE] > *ymax_)
        continue;
    }

    if (passArr[iE] > 0.) {
      pass_.emplace_back(xArr[iE], 0.);
      if (hasY)
        pass_.back().second = yArr[iE];
    }
    else {
      fail_.emplace_back(xArr[iE], 0.);
      if (hasY)
        fail_.back().second = yArr[iE];
    }
  }
}

void
EfficiencyFitter::setLikelihood(RooAbsReal const* _likelihood, RooArgList const* _parameters, RooRealVar* _xvar, RooRealVar* _yvar/* = 0*/)
{
  likelihood_ = _likelihood;

  parameters_.clear();

  xvar_ = _xvar;
  yvar_ = _yvar;

  RooArgSet leaves;
  _likelihood->leafNodeServerList(&leaves);

  if (!leaves.contains(*xvar_) || (yvar_ && !leaves.contains(*yvar_))) {
    std::cerr << "Likelihood function does not contain correct observables" << std::endl;
    likelihood_ = 0;
    xvar_ = 0;
    yvar_ = 0;

    return;
  }

  RooFIter pitr(_parameters->fwdIterator());
  RooAbsArg* param(0);
  while((param = pitr.next())){
    RooRealVar* var(dynamic_cast<RooRealVar*>(param));

    if (!var || !leaves.contains(*var) || var->isConstant())
      continue;
    if (std::find(parameters_.begin(), parameters_.end(), var) == parameters_.end())
      parameters_.push_back(var);
  }
}

int
EfficiencyFitter::fit(int _verbosity/* = 0*/)
{
  if (!likelihood_)
    return -1;

  TVirtualFitter* fitter(TVirtualFitter::Fitter(0, parameters_.size()));
  fitter->Clear();
  fitter->SetFCN(EfficiencyFitter::FCN);

  double verb(_verbosity);
  fitter->ExecuteCommand("SET PRINT", &verb, 1);
  if(_verbosity < 0) fitter->ExecuteCommand("SET NOW", 0, 0);
  else fitter->ExecuteCommand("SET WAR", 0, 0);

  for(unsigned iF(0); iF != parameters_.size(); ++iF) {
    auto* param(parameters_[iF]);
    fitter->SetParameter(iF, param->GetName(), param->getVal(), (param->getMax() - param->getMin()) / 1000., param->getMin(), param->getMax());
  }

  double errdef(1.);
  fitter->ExecuteCommand("SET ERRDEF", &errdef, 1);

  int status(fitter->ExecuteCommand("MINIMIZE", 0, 0));

  if(hesse){
    double maxIter(500. * parameters_.size());
    fitter->ExecuteCommand("HESSE", &maxIter, 1);
  }

  if(minos){
    double maxIter(500. * parameters_.size());
    fitter->ExecuteCommand("MINOS", &maxIter, 1);
  }

  for(unsigned iF(0); iF != parameters_.size(); ++iF){
    char name[100];
    double value;
    double error;
    double vlow;
    double vhigh;
    fitter->GetParameter(iF, name, value, error, vlow, vhigh);
    parameters_[iF]->setVal(value);
    parameters_[iF]->setError(error);
  }

  return status;
}

RooRealVar*
EfficiencyFitter::getParam(char const* _name) const
{
  TString name(_name);
  for(unsigned iF(0); iF != parameters_.size(); ++iF)
    if(parameters_[iF]->GetName() == name) return parameters_[iF];

  return 0;
}

double
EfficiencyFitter::getNLL() const
{
  std::vector<double> parameters(parameters_.size());
  for(unsigned iF(0); iF != parameters_.size(); ++iF)
    parameters[iF] = parameters_[iF]->getVal();

  return fcn(&(parameters[0]));
}

void
EfficiencyFitter::reset()
{
  parameters_.clear();
  likelihood_ = 0;

  xvar_ = 0;
  yvar_ = 0;

  pass_.clear();
  fail_.clear();

  delete xmin_;
  xmin_ = 0;
  delete xmax_;
  xmax_ = 0;
  delete ymin_;
  ymin_ = 0;
  delete ymax_;
  ymax_ = 0;
}

void
EfficiencyFitter::setXRange(double _min, double _max)
{
  if (!xmin_)
    xmin_ = new double(_min);
  else
    *xmin_ = _min;

  if (!xmax_)
    xmax_ = new double(_max);
  else
    *xmax_ = _max;
}

void
EfficiencyFitter::setYRange(double _min, double _max)
{
  if (!ymin_)
    ymin_ = new double(_min);
  else
    *ymin_ = _min;

  if (!ymax_)
    ymax_ = new double(_max);
  else
    *ymax_ = _max;
}

/*static*/
EfficiencyFitter*
EfficiencyFitter::singleton()
{
  static EfficiencyFitter fitter;
  return &fitter;
}
 
/*static*/
void
EfficiencyFitter::FCN(int&, double*, double& _fval, double* _xval, int)
{
  _fval = singleton()->fcn(_xval);
}

double
EfficiencyFitter::fcn(double const* _xval) const
{
  for(unsigned iF(0); iF != parameters_.size(); ++iF)
    parameters_[iF]->setVal(_xval[iF]);

  double nll(0.);

  for (auto& xy : pass_) {
    xvar_->setVal(xy.first);
    if (yvar_)
      yvar_->setVal(xy.second);

    double l(likelihood_->getVal());
    if (l < 0.)
      nll += 1.e+10;
    else
      nll += -std::log(l);
  }

  for (auto& xy : fail_) {
    xvar_->setVal(xy.first);
    if (yvar_)
      yvar_->setVal(xy.second);

    double l(1. - likelihood_->getVal());
    if (l < 0.)
      nll += 1.e+10;
    else
      nll += -std::log(l);
  }

  return nll;
}
