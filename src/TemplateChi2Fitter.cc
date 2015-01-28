#include "../include/TemplateChi2Fitter.h"

#include "TVirtualFitter.h"
#include "TString.h"
#include "TH1D.h"

#include "RooNumber.h"
#include "RooArgList.h"

#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cmath>
#include <set>

TemplateChi2Fitter::TemplateChi2Fitter() :
  nX_(),
  nY_(),
  targetContents_(),
  targetErrors_(),
  templateContents_(),
  templateErrors_(),
  binOffsets_(),
  bins_(),
  chi2_("chi2", "chi2", 0., 0., RooNumber::infinity()),
  scales_(),
  floats_()
{
}

TemplateChi2Fitter::~TemplateChi2Fitter()
{
  for(unsigned iS(0); iS != scales_.size(); ++iS)
    delete scales_[iS];
}

void
TemplateChi2Fitter::setTarget(TH1 const* _hist)
{
  if(!_hist) return;

  nX_.clear();
  nY_.clear();
  targetContents_.clear();
  targetErrors_.clear();
  templateContents_.clear();
  templateErrors_.clear();
  binOffsets_.clear();
  bins_.clear();
  for(unsigned iS(0); iS != scales_.size(); ++iS)
    delete scales_[iS];
  scales_.clear();
  floats_.clear();

  binOffsets_.push_back(0);

  nX_.push_back(_hist->GetNbinsX());
  nY_.push_back(_hist->GetNbinsY());

  for(int iY(1); iY <= nY_[0]; ++iY){
    for(int iX(1); iX <= nX_[0]; ++iX){
      targetContents_.push_back(_hist->GetBinContent(iX, iY));
      targetErrors_.push_back(_hist->GetBinError(iX, iY));

      bins_.push_back(bins_.size());
    }
  }
}

void
TemplateChi2Fitter::setTarget(TObjArray const* _hists)
{
  if(!_hists) return;

  nX_.clear();
  nY_.clear();
  targetContents_.clear();
  targetErrors_.clear();
  templateContents_.clear();
  templateErrors_.clear();
  binOffsets_.clear();
  bins_.clear();
  for(unsigned iS(0); iS != scales_.size(); ++iS)
    delete scales_[iS];
  scales_.clear();
  floats_.clear();

  for(int iH(0); iH != _hists->GetEntries(); ++iH){
    TH1 const* hist(dynamic_cast<TH1 const*>(_hists->At(iH)));
    if(!hist) continue;

    int nX(hist->GetNbinsX());
    int nY(hist->GetNbinsY());

    nX_.push_back(nX);
    nY_.push_back(nY);

    binOffsets_.push_back(bins_.size());

    for(int iY(1); iY <= nY; ++iY){
      for(int iX(1); iX <= nX; ++iX){
        targetContents_.push_back(hist->GetBinContent(iX, iY));
        targetErrors_.push_back(hist->GetBinError(iX, iY));

        bins_.push_back(bins_.size());
      }
    }
  }
}

void
TemplateChi2Fitter::setTarget(int _nB, double const* _cont, double const* _err)
{
  nX_.clear();
  nY_.clear();
  targetContents_.clear();
  targetErrors_.clear();
  templateContents_.clear();
  templateErrors_.clear();
  binOffsets_.clear();
  bins_.clear();
  for(unsigned iS(0); iS != scales_.size(); ++iS)
    delete scales_[iS];
  scales_.clear();
  floats_.clear();

  binOffsets_.push_back(0);

  nX_.push_back(_nB);
  nY_.push_back(1);

  for(int iX(0); iX < nX_[0]; ++iX){
    targetContents_.push_back(_cont[iX]);
    targetErrors_.push_back(_err[iX]);

    bins_.push_back(bins_.size());
  }
}

void
TemplateChi2Fitter::addTemplate(TH1 const* _hist, char const* _name/* = 0*/, RooAbsReal* _scale/* = 0*/)
{
  if(!_hist) return;

  if(_hist->GetNbinsX() != nX_[0] || _hist->GetNbinsY() != nY_[0])
    throw std::runtime_error("Target and template inconsistent");

  TString name(_name);
  if(_scale){
    if(name == "") name = _scale->GetName();
    scales_.push_back(new RooRealProxy(name, TString::Format("scale for %s", _hist->GetName()), &chi2_, *_scale));
    
    RooArgList list;
    _scale->leafNodeServerList(&list);
    for(int iS(0); iS != list.getSize(); ++iS){
      RooRealVar* var(dynamic_cast<RooRealVar*>(list.at(iS)));
      if(!var) continue;
      if(var->isConstant()) continue;
      if(std::find(floats_.begin(), floats_.end(), var) == floats_.end()) floats_.push_back(var);
    }
  }
  else{
    if(name == "") name = _hist->GetName();
    RooRealVar* scale(new RooRealVar(name, name, 1., -RooNumber::infinity(), RooNumber::infinity()));
    scales_.push_back(new RooRealProxy(name, TString::Format("scale for %s", _hist->GetName()), &chi2_, *scale, true, false, true));
    floats_.push_back(scale);

    scales_.back()->Print();
    scales_.back()->arg().Print();
  }

  templateContents_.resize(scales_.size());
  templateErrors_.resize(scales_.size());

  std::vector<double>& contents(templateContents_.back());
  std::vector<double>& errors(templateErrors_.back());

  for(int iY(1); iY <= nY_[0]; ++iY){
    for(int iX(1); iX <= nX_[0]; ++iX){
      contents.push_back(_hist->GetBinContent(iX, iY));
      errors.push_back(_hist->GetBinError(iX, iY));
    }
  }
}

void
TemplateChi2Fitter::addTemplate(TObjArray const* _hists, char const* _name/* = 0*/, RooAbsReal* _scale/* = 0*/)
{
  if(!_hists || _hists->GetEntries() == 0) return;

  TString name(_name);
  if(_scale){
    if(name == "") name = _scale->GetName();
    scales_.push_back(new RooRealProxy(name, TString::Format("scale for %s", _hists->At(0)->GetName()), &chi2_, *_scale));

    RooArgList list;
    _scale->leafNodeServerList(&list);
    for(int iS(0); iS != list.getSize(); ++iS){
      RooRealVar* var(dynamic_cast<RooRealVar*>(list.at(iS)));
      if(!var) continue;
      if(var->isConstant()) continue;
      if(std::find(floats_.begin(), floats_.end(), var) == floats_.end()) floats_.push_back(var);
    }
  }
  else{
    if(name == "") name = _hists->At(0)->GetName();
    RooRealVar* scale(new RooRealVar(name, name, 1. -RooNumber::infinity(), RooNumber::infinity()));
    scales_.push_back(new RooRealProxy(name, "scale for " + name, &chi2_, *scale, true, false, true));
    floats_.push_back(scale);
  }

  templateContents_.resize(scales_.size());
  templateErrors_.resize(scales_.size());

  std::vector<double>& contents(templateContents_.back());
  std::vector<double>& errors(templateErrors_.back());

  unsigned iP(0);

  for(int iH(0); iH != _hists->GetEntries(); ++iH){
    TH1 const* hist(dynamic_cast<TH1 const*>(_hists->At(iH)));
    if(!hist) continue;

    if(hist->GetNbinsX() != nX_[iP] || hist->GetNbinsY() != nY_[iP])
      throw std::runtime_error("Target and template inconsistent");

    for(int iY(1); iY <= nY_[iP]; ++iY){
      for(int iX(1); iX <= nX_[iP]; ++iX){
        contents.push_back(hist->GetBinContent(iX, iY));
        errors.push_back(hist->GetBinError(iX, iY));
      }
    }
    
    ++iP;
  }
}

void
TemplateChi2Fitter::addTemplate(int _nB, double const* _cont, double const* _err, char const* _name/* = 0*/, RooAbsReal* _scale/* = 0*/)
{
  if(_nB != nX_[0])
    throw std::runtime_error("Target and template inconsistent");

  TString name(_name);
  if(_scale){
    if(name == "") name = _scale->GetName();
    scales_.push_back(new RooRealProxy(name, "scale", &chi2_, *_scale));
    
    RooArgList list;
    _scale->leafNodeServerList(&list);
    for(int iS(0); iS != list.getSize(); ++iS){
      RooRealVar* var(dynamic_cast<RooRealVar*>(list.at(iS)));
      if(!var) continue;
      if(var->isConstant()) continue;
      if(std::find(floats_.begin(), floats_.end(), var) == floats_.end()) floats_.push_back(var);
    }
  }
  else{
    if(name == "")
      throw std::runtime_error("No name given to template");
    RooRealVar* scale(new RooRealVar(name, name, 1., -RooNumber::infinity(), RooNumber::infinity()));
    scales_.push_back(new RooRealProxy(name, "scale", &chi2_, *scale, true, false, true));
    floats_.push_back(scale);

    scales_.back()->Print();
    scales_.back()->arg().Print();
  }

  templateContents_.resize(scales_.size());
  templateErrors_.resize(scales_.size());

  std::vector<double>& contents(templateContents_.back());
  std::vector<double>& errors(templateErrors_.back());

  for(int iX(0); iX < nX_[0]; ++iX){
    contents.push_back(_cont[iX]);
    errors.push_back(_err[iX]);
  }
}

void
TemplateChi2Fitter::changeTarget(TH1 const* _hist, unsigned _iH/* = 0*/)
{
  if(!_hist || _iH >= binOffsets_.size()) return;

  unsigned offset(binOffsets_[_iH]);

  unsigned bin(offset);
  for(int iY(1); iY <= nY_[_iH]; ++iY){
    for(int iX(1); iX <= nX_[_iH]; ++iX){
      targetContents_[bin] = _hist->GetBinContent(iX, iY);
      targetErrors_[bin] = _hist->GetBinError(iX, iY);

      ++bin;
    }
  }
}


int
TemplateChi2Fitter::fit(int _verbosity/* = 0*/)
{
  if(targetContents_.size() == 0 || scales_.size() == 0)
    throw std::runtime_error("Fitter empty");

  TVirtualFitter* fitter(TVirtualFitter::Fitter(0, floats_.size()));
  fitter->Clear();
  fitter->SetFCN(TemplateChi2Fitter::FCN);

  double verb(_verbosity);
  fitter->ExecuteCommand("SET PRINT", &verb, 1);
  if(_verbosity < 0) fitter->ExecuteCommand("SET NOW", 0, 0);
  else fitter->ExecuteCommand("SET WAR", 0, 0);

  for(unsigned iF(0); iF != floats_.size(); ++iF)
    fitter->SetParameter(iF, TString("flt_") + floats_[iF]->GetName(), floats_[iF]->getVal(), 0.01, floats_[iF]->getMin(), floats_[iF]->getMax());

  double errdef(1.);
  fitter->ExecuteCommand("SET ERRDEF", &errdef, 1);

  int status(fitter->ExecuteCommand("MINIMIZE", 0, 0));

  for(unsigned iF(0); iF != floats_.size(); ++iF){
    char name[100];
    double value;
    double error;
    double vlow;
    double vhigh;
    fitter->GetParameter(iF, name, value, error, vlow, vhigh);
    floats_[iF]->setVal(value);
    floats_[iF]->setError(error);
  }

  fitter->Clear();
  fitter->SetFCN((void*)NULL);

  return status;
}

RooRealVar*
TemplateChi2Fitter::getFloat(char const* _name) const
{
  TString name(_name);
  for(unsigned iF(0); iF != floats_.size(); ++iF)
    if(floats_[iF]->GetName() == name) return floats_[iF];

  return 0;
}

double
TemplateChi2Fitter::getChi2() const
{
  std::vector<double> floats(floats_.size());
  for(unsigned iF(0); iF != floats_.size(); ++iF)
    floats[iF] = floats_[iF]->getVal();

  return fcn(&(floats[0]));
}

unsigned
TemplateChi2Fitter::getNdof() const
{
  return bins_.size() - floats_.size();
}

void
TemplateChi2Fitter::plot(TDirectory* _dir) const
{
  _dir->cd();

  unsigned nTotalBins(0);
  for(unsigned iP(0); iP != binOffsets_.size(); ++iP)
    nTotalBins += nX_[iP] * nY_[iP];

  TH1D* hTarget(new TH1D("target", "Target", nTotalBins, 0., nTotalBins));
  std::vector<TH1D*> hTemplates;
  for(unsigned iT(0); iT != scales_.size(); ++iT)
    hTemplates.push_back(new TH1D(scales_[iT]->GetName(), scales_[iT]->GetName(), nTotalBins, 0., nTotalBins));

  unsigned iB(0);
  for(unsigned bin(0); bin != nTotalBins; ++bin){
    if(bin != bins_[iB]) continue;

    hTarget->SetBinContent(bin + 1, targetContents_[bin]);
    hTarget->SetBinError(bin + 1, targetErrors_[bin]);

    for(unsigned iT(0); iT != scales_.size(); ++iT){
      hTemplates[iT]->SetBinContent(bin + 1, templateContents_[iT][bin] * (*scales_[iT]));
      hTemplates[iT]->SetBinError(bin + 1, templateErrors_[iT][bin] * (*scales_[iT]));
    }
    ++iB;
  }
}

void
TemplateChi2Fitter::excludeBin(int _iX, int _iY/* = 1*/, int _iP/* = 0*/)
{
  if(nX_[_iP] == 0 || nY_[_iP] == 0) return;

  unsigned bin(-1);
  if(_iY > 1)
    bin = (_iX - 1) + (_iY - 1) * nX_[_iP];
  else
    bin = ((_iX - 1) % nX_[_iP]) + ((_iX - 1) / nX_[_iP]) * nX_[_iP];

  bin += binOffsets_[_iP];

  std::vector<unsigned>::iterator bItr(std::find(bins_.begin(), bins_.end(), bin));
  if(bItr != bins_.end()) bins_.erase(bItr);
}

/*static*/
TemplateChi2Fitter*
TemplateChi2Fitter::singleton()
{
  static TemplateChi2Fitter fitter;
  return &fitter;
}

/*static*/
void
TemplateChi2Fitter::FCN(int&, double*, double& _fval, double* _xval, int)
{
  _fval = singleton()->fcn(_xval);
}

double
TemplateChi2Fitter::fcn(double const* _xval) const
{
  for(unsigned iF(0); iF != floats_.size(); ++iF)
    floats_[iF]->setVal(_xval[iF]);
  
  double chi2(0.);

  unsigned nT(scales_.size());
  unsigned nB(bins_.size());

  std::vector<double> scales(nT);
  for(unsigned iT(0); iT != nT; ++iT)
    scales[iT] = scales_[iT]->arg().getVal();

  for(unsigned iB(0); iB != nB; ++iB){
    unsigned bin(bins_[iB]);

    double rtnumer(targetContents_[bin]);
    double denom(targetErrors_[bin] * targetErrors_[bin]);
    for(unsigned iT(0); iT != nT; ++iT){
      rtnumer -= templateContents_[iT][bin] * scales[iT];
      denom += templateErrors_[iT][bin] * templateErrors_[iT][bin] * scales[iT] * scales[iT];
    }
    if(denom == 0.) continue;

    chi2 += rtnumer * rtnumer / denom;
  }

  return chi2;
}
