#include "../include/IntegralChi2Fitter.h"

#include "TVirtualFitter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "RooAbsPdf.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstring>

IntegralChi2Fitter::IntegralChi2Fitter() :
  hesse(false),
  minos(false),
  dimension_(0),
  functionName_(""),
  binIndices_(),
  bins_(),
  covariance_(0),
  observables_(),
  extendPdf_(0)
{
}

IntegralChi2Fitter::~IntegralChi2Fitter()
{
  reset();
}

void
IntegralChi2Fitter::setTarget(TH1 const* _hist, bool _normalized/* = true*/, char const* _binNameSuffix/* = 0*/)
{
  if(!_hist) return;

  dimension_ = _hist->GetDimension();

  reset();

  Bin bin(dimension_);

  binEdges_.assign(dimension_, std::vector<double>());
  binEdges_[0].resize(_hist->GetNbinsX() + 1);
  if(dimension_ > 1) binEdges_[1].resize(_hist->GetNbinsY() + 1);
  if(dimension_ > 2) binEdges_[2].resize(_hist->GetNbinsZ() + 1);

  for(int iZ(1); iZ <= _hist->GetNbinsZ(); ++iZ){
    if(dimension_ > 2){
      binEdges_[2][iZ - 1] = _hist->GetZaxis()->GetBinLowEdge(iZ);
      binEdges_[2][iZ] = _hist->GetZaxis()->GetBinUpEdge(iZ);
      bin.setIndex(2, iZ - 1);
    }

    for(int iY(1); iY <= _hist->GetNbinsY(); ++iY){
      if(dimension_ > 1){
        binEdges_[1][iY - 1] = _hist->GetYaxis()->GetBinLowEdge(iY);
        binEdges_[1][iY] = _hist->GetYaxis()->GetBinUpEdge(iY);
        bin.setIndex(1, iY - 1);
      }

      for(int iX(1); iX <= _hist->GetNbinsX(); ++iX){
        binEdges_[0][iX - 1] = _hist->GetXaxis()->GetBinLowEdge(iX);
        binEdges_[0][iX] = _hist->GetXaxis()->GetBinUpEdge(iX);
        bin.setIndex(0, iX - 1);

        int binIndex(binIndices_.size());
        binIndices_.push_back(binIndex);

        bin.name = TString::Format("bin%d", binIndex);
        if(_binNameSuffix) bin.name += _binNameSuffix;

        bin.volume = binEdges_[0][iX] - binEdges_[0][iX - 1];
        if(dimension_ > 1) bin.volume *= binEdges_[1][iY] - binEdges_[1][iY - 1];
        if(dimension_ > 2) bin.volume *= binEdges_[2][iZ] - binEdges_[2][iZ - 1];

        double norm(_normalized ? bin.volume : 1.);

        int iBin(_hist->GetBin(iX, iY, iZ));

        bin.content = _hist->GetBinContent(iBin) * norm;
        bin.error = _hist->GetBinError(iBin) * norm;

        bins_.push_back(bin);
      }
    }
  }
}

void
IntegralChi2Fitter::setTarget(int _nB, double const* _edges, double const* _cont, double const* _err, bool _normalized/* = true*/, char const* _binNameSuffix/* = 0*/)
{
  dimension_ = 1;

  reset();

  Bin bin(dimension_);

  binEdges_.assign(1, std::vector<double>(_edges, _edges + _nB + 1));

  for(int iX(0); iX < _nB; ++iX){
    if(_edges[iX + 1] <= _edges[iX])
      throw std::runtime_error("Bin edges not in increasing order");

    bin.setIndex(0, iX);

    binIndices_.push_back(iX);

    bin.name = TString::Format("bin%d", iX);
    if(_binNameSuffix) bin.name += _binNameSuffix;

    bin.volume = binEdges_[0][iX + 1] - binEdges_[0][iX];

    double norm(_normalized ? bin.volume : 1.);

    bin.content = _cont[iX] * norm;
    bin.error = _err[iX] * norm;

    bins_.push_back(bin);
  }
}

void
IntegralChi2Fitter::setTarget(int _nB, double _low, double _high, double const* _cont, double const* _err, bool _normalized/* = true*/, char const* _binNameSuffix/* = 0*/)
{
  dimension_ = 1;

  reset();

  Bin bin(dimension_);

  binEdges_.assign(1, std::vector<double>(_nB + 1));

  double delta((_high - _low) / _nB);

  for(int iX(0); iX < _nB; ++iX){
    binEdges_[0][iX] = _low + delta * iX;
    binEdges_[0][iX + 1] = _low + delta * (iX + 1);

    bin.setIndex(0, iX);

    binIndices_.push_back(iX);

    bin.name = TString::Format("bin%d", iX);
    if(_binNameSuffix) bin.name += _binNameSuffix;

    bin.volume = binEdges_[0][iX + 1] - binEdges_[0][iX];

    double norm(_normalized ? bin.volume : 1.);

    bin.content = _cont[iX] * norm;
    bin.error = _err[iX] * norm;

    bins_.push_back(bin);
  }
}

void
IntegralChi2Fitter::setFunction(RooAbsReal const* _func, RooRealVar* _observable, char const* _sliceName/* = 0*/, RooArgSet* _slicedVars/* = 0*/)
{
  if(!_observable) return;

  RooArgList args(*_observable);
  setFunction(_func, &args, _sliceName, _slicedVars);
}

void
IntegralChi2Fitter::setFunction(RooAbsReal const* _func, RooArgList const* _observables, char const* _sliceName/* = 0*/, RooArgSet* _slicedVars/* = 0*/)
{
  if(_observables->getSize() != dimension_)
    throw std::logic_error("Function and target have different dimensionalities");
  for(int iD(0); iD != dimension_; ++iD){
    if(!_observables->at(iD)->InheritsFrom(RooRealVar::Class()))
      throw std::runtime_error("Observable is not a RooRealVar");
  }

  if(binIndices_.size() == 0 || !_func || !_observables) return;

  clearFunction();

  functionName_ = _func->GetName();

  observables_.add(*_observables);

  if(_func->InheritsFrom(RooAbsPdf::Class())){
    extendPdf_ = static_cast<RooAbsPdf const*>(_func);
    if(!extendPdf_->canBeExtended())
      throw std::runtime_error("Non-extendable PDF passed as fit function");

    std::cout << "Using extended PDF as fit function" << std::endl;
  }

  RooArgSet leaves;
  _func->leafNodeServerList(&leaves);
  RooFIter litr(leaves.fwdIterator());
  RooAbsArg* leaf(0);
  while((leaf = litr.next())){
    RooRealVar* var(dynamic_cast<RooRealVar*>(leaf));
    // contains only checks the identity of the names
    if(!var || _observables->contains(*var) || var->isConstant()) continue;
    if(std::find(params_.begin(), params_.end(), var) == params_.end()) params_.push_back(var);
  }

  for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex){
    Bin& bin(bins_[binIndex]);
    for(int iD(0); iD != dimension_; ++iD){
      unsigned index(bin.index(iD));
      static_cast<RooRealVar*>(_observables->at(iD))->setRange(bin.name, binEdges_[iD][index], binEdges_[iD][index + 1]);
    }

    if(_sliceName && _slicedVars){
      RooFIter litr(_slicedVars->fwdIterator());
      RooAbsArg* leaf(0);
      while((leaf = litr.next())){
        RooRealVar* var(dynamic_cast<RooRealVar*>(leaf));
        if(var && var->hasBinning(_sliceName)){
          double min(var->getMin(_sliceName));
          double max(var->getMax(_sliceName));
          var->setRange(bin.name, min, max);
        }
      }
    }

    // If slice / projection is specified, need to integrated over those variables too.
    RooArgSet iset(*_observables);
    if(_slicedVars) iset.add(*_slicedVars);

    if(extendPdf_){
      // Setting the observables as "normalization set" - component pdfs of the func will be normalized to 1 over the full range of the observables.
      // This option is not completely thought through - might break something..
      bin.integral = _func->createIntegral(iset, iset, bin.name);
    }
    else
      bin.integral = _func->createIntegral(iset, bin.name);

    bin.integral->SetName(bin.name + "_integral");
  }
}

int
IntegralChi2Fitter::fit(int _verbosity/* = 0*/)
{
  TVirtualFitter* fitter(TVirtualFitter::Fitter(0, params_.size()));
  fitter->Clear();
  fitter->SetFCN(IntegralChi2Fitter::FCN);

  double verb(_verbosity);
  fitter->ExecuteCommand("SET PRINT", &verb, 1);
  if(_verbosity < 0) fitter->ExecuteCommand("SET NOW", 0, 0);
  else fitter->ExecuteCommand("SET WAR", 0, 0);

  for(unsigned iF(0); iF != params_.size(); ++iF)
    fitter->SetParameter(iF, params_[iF]->GetName(), params_[iF]->getVal(), 0.01, 0., 0.);

  double errdef(1.);
  fitter->ExecuteCommand("SET ERRDEF", &errdef, 1);

  int status(fitter->ExecuteCommand("MINIMIZE", 0, 0));

  if(hesse){
    double maxIter(500. * params_.size());
    fitter->ExecuteCommand("HESSE", &maxIter, 1);
  }

  if(minos){
    double maxIter(500. * params_.size());
    fitter->ExecuteCommand("MINOS", &maxIter, 1);
  }

  for(unsigned iF(0); iF != params_.size(); ++iF){
    char name[100];
    double value;
    double error;
    double vlow;
    double vhigh;
    fitter->GetParameter(iF, name, value, error, vlow, vhigh);
    params_[iF]->setVal(value);
    params_[iF]->setError(error);
  }

  covariance_ = new TMatrixDSym(params_.size(), fitter->GetCovarianceMatrix());

  return status;
}

RooRealVar*
IntegralChi2Fitter::getParam(char const* _name) const
{
  TString name(_name);
  for(unsigned iF(0); iF != params_.size(); ++iF)
    if(params_[iF]->GetName() == name) return params_[iF];

  return 0;
}

double
IntegralChi2Fitter::getChi2() const
{
  std::vector<double> params(params_.size());
  for(unsigned iF(0); iF != params_.size(); ++iF)
    params[iF] = params_[iF]->getVal();

  return fcn(&(params[0]));
}

unsigned
IntegralChi2Fitter::getNdof() const
{
  return binIndices_.size() - params_.size();
}

void
IntegralChi2Fitter::excludeBin(int _iX)
{
  std::vector<int>::iterator itr(std::find(binIndices_.begin(), binIndices_.end(), _iX));
  if(itr != binIndices_.end())
    binIndices_.erase(itr);
}

void
IntegralChi2Fitter::reset()
{
  clearFunction();
  binIndices_.clear();
  bins_.clear();
  binEdges_.clear();
}

TH1*
IntegralChi2Fitter::histogram() const
{
  if(bins_.size() == 0) return 0;
  if(dimension_ > 3) return 0;

  TH1* histo(0);

  if(dimension_ == 1)
    histo = new TH1D(functionName_ + "_integral", "Integral of " + functionName_, binEdges_[0].size() - 1, &binEdges_[0][0]);
  else if(dimension_ == 2)
    histo = new TH2D(functionName_ + "_integral", "Integral of " + functionName_, binEdges_[0].size() - 1, &binEdges_[0][0], binEdges_[1].size() - 1, &binEdges_[1][0]);
  else if(dimension_ == 3)
    histo = new TH3D(functionName_ + "_integral", "Integral of " + functionName_, binEdges_[0].size() - 1, &binEdges_[0][0], binEdges_[1].size() - 1, &binEdges_[1][0], binEdges_[2].size() - 1, &binEdges_[2][0]);

  for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex){
    Bin const& bin(bins_[binIndex]);
    int iBin(histo->GetBin(bin.index(0) + 1, bin.index(1) + 1, bin.index(2) + 1));

    double content(bin.integral->getVal() / bin.volume);
    if(extendPdf_) content *= extendPdf_->expectedEvents(observables_);

    histo->SetBinContent(iBin, content);
  }

  return histo;
}

TH1*
IntegralChi2Fitter::pullHistogram() const
{
  TH1* histo(histogram());
  if(!histo) return 0;

  histo->SetTitle("Pull (target - " + functionName_ + ")");

  TH1* target(static_cast<TH1*>(histo->Clone("target")));
  TH1* errors(static_cast<TH1*>(histo->Clone("errors")));

  for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex){
    Bin const& bin(bins_[binIndex]);
    int iBin(histo->GetBin(bin.index(0) + 1, bin.index(1) + 1, bin.index(2) + 1));

    target->SetBinContent(iBin, bin.content / bin.volume);
    errors->SetBinContent(iBin, bin.error / bin.volume);
  }

  histo->Add(target, -1.);
  histo->Scale(-1.);
  histo->Divide(errors);

  delete target;
  delete errors;

  return histo;
}

double
IntegralChi2Fitter::integral() const
{
  double sum(0.);
  for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex)
    if(bins_[binIndex].integral) sum += bins_[binIndex].integral->getVal();

  if(extendPdf_) sum *= extendPdf_->expectedEvents(observables_);

  return sum;
}

double
IntegralChi2Fitter::sumTarget() const
{
  double sum(0.);
  for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex)
    sum += bins_[binIndex].content;

  return sum;
}

/*static*/
IntegralChi2Fitter*
IntegralChi2Fitter::singleton()
{
  static IntegralChi2Fitter fitter;
  return &fitter;
}
 
/*static*/
void
IntegralChi2Fitter::FCN(int&, double*, double& _fval, double* _xval, int)
{
  _fval = singleton()->fcn(_xval);
}

void
IntegralChi2Fitter::clearFunction()
{
  // List of binning name is a global. Addtionally removeRange only sets the range limits to -inf, inf and not delete it
  RooFIter oitr(observables_.fwdIterator());
  RooRealVar* var(0);
  while((var = static_cast<RooRealVar*>(oitr.next()))){
    double defaultMax(var->getMax());
    double defaultMin(var->getMin());
    for(unsigned binIndex(0); binIndex != bins_.size(); ++binIndex)
      var->setRange(bins_[binIndex].name, defaultMin, defaultMax);
  }

  extendPdf_ = 0;
  observables_.removeAll();
  functionName_ = "";

  params_.clear();
  delete covariance_;
  covariance_ = 0;
}

double
IntegralChi2Fitter::fcn(double const* _xval) const
{
  double chi2(0.);

  for(unsigned iF(0); iF != params_.size(); ++iF)
    params_[iF]->setVal(_xval[iF]);

  double scale(extendPdf_ ? extendPdf_->expectedEvents(observables_) : 1.);

  for(unsigned iB(0); iB != binIndices_.size(); ++iB){
    Bin const& bin(bins_[binIndices_[iB]]);
    double denom(bin.error * bin.error);
    if(denom == 0.) continue;

    double rtnumer(bin.content - bin.integral->getVal() * scale);

    chi2 += rtnumer * rtnumer / denom;
  }

  return chi2;
}


IntegralChi2Fitter::Bin::Bin() :
  ndim(0),
  content(0.),
  error(0.),
  volume(0.),
  integral(0),
  index_(0)
{
}

IntegralChi2Fitter::Bin::Bin(unsigned _ndim) :
  ndim(_ndim),
  name(""),
  content(0.),
  error(0.),
  volume(0.),
  integral(0),
  index_(new unsigned[_ndim])
{
}

IntegralChi2Fitter::Bin::Bin(IntegralChi2Fitter::Bin const& _orig) :
  ndim(_orig.ndim),
  name(_orig.name),
  content(_orig.content),
  error(_orig.error),
  volume(_orig.volume),
  integral(_orig.integral ? static_cast<RooAbsReal*>(_orig.integral->Clone()) : 0),
  index_(new unsigned[_orig.ndim])
{
  for(unsigned iD(0); iD != ndim; ++iD)
    index_[iD] = _orig.index_[iD];
}
  
IntegralChi2Fitter::Bin::~Bin()
{
  delete integral;
  delete [] index_;
}

IntegralChi2Fitter::Bin&
IntegralChi2Fitter::Bin::operator=(IntegralChi2Fitter::Bin const& _rhs)
{
  delete integral;
  delete [] index_;

  ndim = _rhs.ndim;
  name = _rhs.name;
  content = _rhs.content;
  error = _rhs.error;
  volume = _rhs.volume;
  integral = _rhs.integral ? static_cast<RooAbsReal*>(_rhs.integral->Clone()) : 0;
  index_ = new unsigned[ndim];

  return *this;
}

void
IntegralChi2Fitter::Bin::setIndex(unsigned _iD, int _iX)
{
  if(_iD >= ndim) return;
  index_[_iD] = _iX;
}

int
IntegralChi2Fitter::Bin::index(unsigned _iD) const
{
  if(_iD >= ndim) return 0;
  else return index_[_iD];
}
