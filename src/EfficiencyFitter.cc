#include "../include/EfficiencyFitter.h"

#include "TTree.h"
#include "TVirtualFitter.h"
#include "RooArgSet.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

EfficiencyFitter::~EfficiencyFitter()
{
  reset();
}

void
EfficiencyFitter::setTarget(TTree* _tree, char const* _passdef, char const* _xname, char const* _yname/* = 0*/)
{
  if(!_tree) return;

  reset();

  xname_ = _xname;
  yname_ = _yname;

  bool hasY(yname_.Length() != 0);
  
  TString expr(_passdef);
  expr += " != 0:" + xname_;
  if (hasY)
    expr += ":" + yname_;

  _tree->SetEstimate(_tree->GetEntries() + 1);
  long nEntries(_tree->Draw(expr, "", "goff"));
  double* passArr(_tree->GetV1());
  double* xArr(_tree->GetV2());
  double* yArr(_tree->GetV3());

  for (long iE(0); iE != nEntries; ++iE) {
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
EfficiencyFitter::setLikelihood(RooAbsReal const* _likelihood, RooArgList const* _parameters)
{
  likelihood_ = _likelihood;

  parameters_.clear();

  RooArgSet leaves;
  _likelihood->leafNodeServerList(&leaves);

  xvar_ = dynamic_cast<RooRealVar*>(leaves.find(xname_));
  if (yname_.Length() != 0)
    yvar_ = dynamic_cast<RooRealVar*>(leaves.find(yname_));

  RooFIter pitr(_parameters->fwdIterator());
  RooAbsArg* param(0);
  while((param = pitr.next())){
    RooRealVar* var(dynamic_cast<RooRealVar*>(param));

    if (!var || !leaves.contains(*var) || var->isConstant())
      continue;
    if (std::find(parameters_.begin(), parameters_.end(), var) == parameters_.end())
      parameters_.push_back(var);
  }

  if (!xvar_ || (yname_.Length() != 0 && !yvar_)) {
    std::cerr << "Likelihood function does not contain correct observables" << std::endl;
    likelihood_ = 0;
    xvar_ = 0;
    yvar_ = 0;
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

  for(unsigned iF(0); iF != parameters_.size(); ++iF)
    fitter->SetParameter(iF, parameters_[iF]->GetName(), parameters_[iF]->getVal(), 0.01, 0., 0.);

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

  xname_ = "";
  yname_ = "";
  xvar_ = 0;
  yvar_ = 0;

  pass_.clear();
  fail_.clear();
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

    nll += -std::log(likelihood_->getVal());
  }

  for (auto& xy : fail_) {
    xvar_->setVal(xy.first);
    if (yvar_)
      yvar_->setVal(xy.second);

    nll += -std::log(1. - likelihood_->getVal());
  }

  return nll;
}
