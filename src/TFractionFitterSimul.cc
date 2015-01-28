#include "../include/TFractionFitterSimul.h"

#include "Riostream.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TFitter.h"
#include "TMath.h"
#include "TH2.h"

#include <stdexcept>
#include <iostream>

extern TVirtualFitter* fractionFitter;

TObjArray TFractionFitterSimul::fgFitters;

ClassImp(TFractionFitterExt)
ClassImp(TFractionFitterSimul)

TFractionFitterExt::TFractionFitterExt() :
  TFractionFitter()
{
}

TFractionFitterExt::TFractionFitterExt(TH1* _data, TObjArray* _MCs, Option_t* _opt/* = ""*/) :
  TFractionFitter(_data, _MCs, _opt)
{
  for(int ix(1); ix <= _data->GetNbinsX(); ++ix){
    for(int iy(1); iy <= _data->GetNbinsY(); ++iy){
      for(int iz(1); iz <= _data->GetNbinsZ(); ++iz){
        int bin(_data->GetBin(ix, iy, iz));
	unsigned nEmpty(0);
        for(int iP(0); iP < _MCs->GetEntries(); ++iP)
          if(static_cast<TH1*>(_MCs->At(iP))->GetBinContent(bin) == 0) ++nEmpty;

        if(nEmpty > 1) ExcludeBin(bin);
      }
    }
  }
}

TFractionFitterExt::~TFractionFitterExt()
{
}

// complete copy of TFractionFitter::ComputeChisquareLambda (private)
void
TFractionFitterExt::ComputeChisquareLambda()
{
  if ( !fFitDone ) {
    Error("ComputeChisquareLambda","Fit not yet (successfully) performed");
    fChisquare = 0;
    return;
  }

  Int_t minX, maxX, minY, maxY, minZ, maxZ;
  GetRanges(minX, maxX, minY, maxY, minZ, maxZ);

  Double_t logLyn = 0; // likelihood of prediction
  Double_t logLmn = 0; // likelihood of data ("true" distribution)
  for(Int_t x = minX; x <= maxX; x++) {
    for(Int_t y = minY; y <= maxY; y++) {
      for(Int_t z = minZ; z <= maxZ; z++) {
        if (IsExcluded(fData->GetBin(x, y, z))) continue;
        Double_t di = fData->GetBinContent(x, y, z);
        Double_t fi = fPlot->GetBinContent(x, y, z);
        if(fi != 0) logLyn += di * TMath::Log(fi) - fi;
        if(di != 0) logLmn += di * TMath::Log(di) - di;
        for(Int_t j = 0; j < fNpar; j++) {
          Double_t aji = ((TH1*)fMCs.At(j))->GetBinContent(x, y, z);
          Double_t bji = ((TH1*)fAji.At(j))->GetBinContent(x, y, z);
          if(bji != 0) logLyn += aji * TMath::Log(bji) - bji;
          if(aji != 0) logLmn += aji * TMath::Log(aji) - aji;
        }
      }
    }
  }

  fChisquare = -2*logLyn + 2*logLmn;
}

// complete copy of TFractionFitter::IsExcluded (private)
Bool_t
TFractionFitterExt::IsExcluded(Int_t bin) const
{
  // Function for internal use, checking whether the given bin is
  // excluded from the fit or not.

  for (unsigned int b = 0; b < fExcludedBins.size(); ++b) 
    if (fExcludedBins[b] == bin) return kTRUE;
  return kFALSE;
}

// complete copy of TFractionFitter::GetRanges (private)
void
TFractionFitterExt::GetRanges(Int_t& minX, Int_t& maxX, Int_t& minY, Int_t& maxY, Int_t& minZ, Int_t& maxZ) const
{
   // Used internally to obtain the bin ranges according to the dimensionality of
   // the histogram and the limits set by hand.

   if (fData->GetDimension() < 2) {
      minY = maxY = minZ = maxZ = 0;
      minX = fLowLimitX;
      maxX = fHighLimitX;
   } else if (fData->GetDimension() < 3) {
      minZ = maxZ = 0;
      minX = fLowLimitX;
      maxX = fHighLimitX;
      minY = fLowLimitY;
      maxY = fHighLimitY;
   } else {
      minX = fLowLimitX;
      maxX = fHighLimitX;
      minY = fLowLimitY;
      maxY = fHighLimitY;
      minZ = fLowLimitZ;
      maxZ = fHighLimitZ;
   }
}

TH1*
TFractionFitterExt::GetPlot() const
{
  if(!fFitDone){
    Error("GetPlot", "Fit not yet performed");
    return 0;
  }

  return fPlot;
}

TFractionFitterSimul::TFractionFitterSimul() :
  TObject(),
  fNpar(0),
  fFitters(),
  fFitterImp(),
  fDirectories()
{
  fFitters.SetOwner(kTRUE);
  fDirectories.SetOwner(kTRUE);
}

TFractionFitterSimul::TFractionFitterSimul(TH1* _data, TObjArray* _MCs, Option_t* _opt/* = ""*/) :
  TObject(),
  fNpar(0),
  fFitters(),
  fFitterImp(),
  fDirectories()
{
  if(!_data || !_MCs || _MCs->GetEntries() == 0){
    Error("TFractionFitterSimul", "Invalid input to constructor");
    return;
  }

  fFitters.SetOwner(kTRUE);
  fDirectories.SetOwner(kTRUE);

  fNpar = _MCs->GetEntries();

  TString opt(_opt);
  opt.ToUpper();

  Bool_t xProj(opt.Contains("PX"));
  Bool_t yProj(opt.Contains("PY"));

  if(xProj || yProj){
    if(_data->GetDimension() == 2){
      TH2* data(static_cast<TH2*>(_data));

      Int_t nSlices(xProj ? data->GetNbinsY() : data->GetNbinsX());

      TDirectory* pwd(gDirectory);

      if(pwd->InheritsFrom(TDirectoryFile::Class())){
        Error("TFractionFitterSimul", "TFractionFitterSimul needs to operate in a non-file directory");
        return;
      }

      for(Int_t iSlice(1); iSlice <= nSlices; ++iSlice){
        TString sVar(TString::Format("var%d", iSlice - 1));
        TDirectory* dir(pwd->mkdir("TFractionFitterSimul_" + sVar, "TFractionFitterSimul " + sVar));
        if(!dir){
          Error("TFractionFitterSimul", "Could not construct directory for " + sVar);
          return;
        }

        fDirectories.Add(dir);

        TString suffix(TString::Format("_%d", iSlice));
        TH1* dataSlice;
        if(xProj) dataSlice = data->ProjectionX(data->GetName() + suffix, iSlice, iSlice);
        else dataSlice = data->ProjectionY(data->GetName() + suffix, iSlice, iSlice);

        dataSlice->SetDirectory(dir);

        TObjArray MCs;
        for(Int_t iPar(0); iPar < fNpar; ++iPar){
          TH2* MC(dynamic_cast<TH2*>(_MCs->At(iPar)));
          if(!MC){
            Error("TFractionFitterSimul", "Invalid input to constructor");
            return;
          }
          TH1* MCSlice;
          if(xProj) MCSlice = MC->ProjectionX(MC->GetName() + suffix, iSlice, iSlice);
          else MCSlice = MC->ProjectionY(MC->GetName() + suffix, iSlice, iSlice);
          MCs.Add(MCSlice);

          MCSlice->SetDirectory(dir);
        }

        AddDataAndMCs(dataSlice, &MCs);
      }
    }

    return;
  }

  AddDataAndMCs(_data, _MCs);
}

TFractionFitterSimul::TFractionFitterSimul(TObjArray const& _dataArray, TObjArray const& _MCsArray, Option_t* _opt/* = ""*/) :
  TObject(),
  fNpar(0),
  fFitters(),
  fFitterImp(),
  fDirectories()
{
  if(_dataArray.GetEntries() == 0 || _MCsArray.GetEntries() == 0 || _dataArray.GetEntries() != _MCsArray.GetEntries()){
    Error("TFractionFitterSimul", "Invalid input to constructor");
    return;
  }
  TObjArray* MCs(dynamic_cast<TObjArray*>(_MCsArray.At(0)));
  if(!MCs){
    Error("TFractionFitterSimul", "Invalid input to constructor");
    return;
  }

  fFitters.SetOwner(kTRUE);
  fDirectories.SetOwner(kTRUE);

  fNpar = MCs->GetEntries();

  TString opt(_opt);
  opt.ToUpper();
  // slicing not implemented

  AddDataAndMCs(_dataArray, _MCsArray);
}

TFractionFitterSimul::~TFractionFitterSimul()
{
  // because TFractionFitter cannot take care of itself
  delete fractionFitter;
  fractionFitter = 0;
}

void
TFractionFitterSimul::AddDataAndMCs(TH1* _data, TObjArray* _MCs)
{
  if(!_data || !_MCs || _MCs->GetEntries() != fNpar){
    Error("AddDataAndMCs", "Invalid input");
    return;
  }

  TDirectory* pwd(gDirectory);

  if(pwd->InheritsFrom(TDirectoryFile::Class())){
    Error("AddDataAndMCs", "TFractionFitterSimul needs to operate in a non-file directory");
    return;
  }

  TString sVar(TString::Format("var%d", fFitters.GetEntries()));
  TString dirName("TFractionFitterSimul_" + sVar);

  TDirectory* dir(static_cast<TDirectory*>(fDirectories.FindObject(dirName)));

  if(!dir){
    dir = pwd->mkdir(dirName, "TFractionFitterSimul " + sVar);
    if(!dir){
      Error("AddDataAndMCs", "Could not construct directory for " + sVar);
      return;
    }
    fDirectories.Add(dir);
  }

  dir->cd();
  TFractionFitterExt* fitter(new TFractionFitterExt(_data, _MCs, "Q"));
  fFitters.Add(fitter);

  pwd->cd();
}

void
TFractionFitterSimul::AddDataAndMCs(TObjArray const& _dataArray, TObjArray const& _MCsArray)
{
  if(_dataArray.GetEntries() == 0 || _MCsArray.GetEntries() == 0 || _dataArray.GetEntries() != _MCsArray.GetEntries()){
    Error("AddDataAndMCs", "Invalid input");
    return;
  }

  for(Int_t iVar(0); iVar != _dataArray.GetEntries(); ++iVar){
    TH1* data(dynamic_cast<TH1*>(_dataArray.At(iVar)));
    TObjArray* MCs(dynamic_cast<TObjArray*>(_MCsArray.At(iVar)));
    if(!data || !MCs || MCs->GetEntries() != fNpar){
      Error("AddDataAndMCs", "Invalid input");
      return;
    }

    AddDataAndMCs(data, MCs);
  }
}

void
TFractionFitterSimul::Constrain(Int_t _iPar, Double_t _low, Double_t _high)
{
  if(_iPar < 0 || _iPar >= fNpar){
    Error("Constrain", "Wrong parameter number");
    return;
  }

  Double_t plist[3] = {Double_t(_iPar), _low, _high};
  fFitterImp.ExecuteCommand("SET LIMIT", plist, 3);
}

void
TFractionFitterSimul::ErrorAnalysis(Double_t _UP)
{
  if(!GetFitDone()){
    Error("ErrorAnalysis","Fit not yet performed");
    return;
  }

  Double_t param(_UP > 0. ? _UP : 0.5);
  fFitterImp.ExecuteCommand("SET ERRDEF", &param, 1);

  Int_t status(fFitterImp.ExecuteCommand("MINOS", 0, 0));
  if(status != 0){
    Error("ErrorAnalysis", "Error return from MINOS: %d",status);
  }
}

Int_t
TFractionFitterSimul::Fit(Option_t* _opt/* = ""*/)
{
  fFitterImp.SetObjectFit(this);
  fFitterImp.SetFCN(TFractionFitterSimul::FCN);

  fgFitters = fFitters;

  // set print level 
  TString opt(_opt);
  opt.ToUpper();
  Double_t plist;
  if (opt.Contains("Q") ) { 
    plist = -1;
    fFitterImp.ExecuteCommand("SET PRINT", &plist, 1);
    fFitterImp.ExecuteCommand("SET NOW", &plist, 0);
  }
  else if (opt.Contains("V") ) { 
    plist = 1;
    fFitterImp.ExecuteCommand("SET PRINT", &plist, 1);
  }

  Double_t defaultFraction = 1.0/((Double_t)fNpar);
  Double_t defaultStep = 0.01;
  for(Int_t par = 0; par < fNpar; ++par) {
    TString name("frac"); name += par;
    fFitterImp.SetParameter(par, name.Data(), defaultFraction, defaultStep, 0, 0);
  }

  plist = 0.5;
  // set the UP value to 0.5
  fFitterImp.ExecuteCommand("SET ERRDEF", &plist, 1);

  // fit
  Int_t status(fFitterImp.ExecuteCommand("MINIMIZE",0,0));

  if (status == 0){
    Int_t dummy1;
    Double_t dummy2;
    Double_t dummy3;
    Double_t dummy4;
    char dummy5[100];
    Double_t* fractions(new Double_t[fNpar]);

    for(Int_t iPar(0); iPar != fNpar; ++iPar)
      fFitterImp.GetParameter(iPar, dummy5, fractions[iPar], dummy2, dummy3, dummy4);

    FCN(dummy1, &dummy2, dummy3, fractions, 3);

    delete fractions;

    for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
      TFractionFitterExt* fitter(static_cast<TFractionFitterExt*>(fFitters.At(iVar)));
      fitter->SetFitDone(kTRUE);
    }

    // determine goodness of fit
    for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
      TFractionFitterExt* fitter(static_cast<TFractionFitterExt*>(fFitters.At(iVar)));
      fitter->ComputeChisquareLambda();
    }
  }

  return status;
}

Double_t
TFractionFitterSimul::GetChisquare() const
{
  Double_t chi2(0.);

  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
    TFractionFitter const* fitter(static_cast<TFractionFitter const*>(fFitters.At(iVar)));
    chi2 += fitter->GetChisquare();
  }

  return chi2;
}

Bool_t
TFractionFitterSimul::GetFitDone() const
{
  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
    TFractionFitterExt const* fitter(static_cast<TFractionFitterExt const*>(fFitters.At(iVar)));
    if(!fitter->GetFitDone()) return kFALSE;
  }

  return kTRUE;
}

TObjArray
TFractionFitterSimul::GetMCPredictions(Int_t _iPar) const
{
  TObjArray predictions;
  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
    TH1* prediction(static_cast<TFractionFitter*>(fFitters.At(iVar))->GetMCPrediction(_iPar));
    if(!prediction){
      predictions.Clear();
      return predictions;
    }
    predictions.Add(prediction);
  }

  return predictions;
}

Int_t
TFractionFitterSimul::GetNDF() const
{
  Int_t result(0);

  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar)
    result += static_cast<TFractionFitter*>(fFitters.At(iVar))->GetNDF();

  return result;
}

TObjArray
TFractionFitterSimul::GetPlots() const
{
  TObjArray plots;

  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar){
    TH1* plot(static_cast<TFractionFitterExt*>(fFitters.At(iVar))->GetPlot());
    if(!plot){
      plots.Clear();
      return plots;
    }
    plots.Add(plot);
  }

  return plots;
}

Double_t
TFractionFitterSimul::GetProb() const
{
  Double_t result(1.);

  for(Int_t iVar(0); iVar != fFitters.GetEntries(); ++iVar)
    result *= static_cast<TFractionFitter*>(fFitters.At(iVar))->GetProb();

  return result;
}

void
TFractionFitterSimul::GetResult(Int_t _iPar, Double_t& _value, Double_t& _error) const
{
  if(_iPar < 0 || _iPar >= fNpar){
    Error("GetResult", "Wrong parameter number");
    return;
  }
  if(!GetFitDone()){
    Error("GetResult","Fit not yet performed");
    return;
  }
  char parname[100];
  Double_t vlow, vhigh;

  fFitterImp.GetParameter(_iPar, parname, _value, _error, vlow, vhigh);
}

void
TFractionFitterSimul::UnConstrain(Int_t _iPar)
{
  if(_iPar < 0 || _iPar >= fNpar){
    Error("Constrain", "Wrong parameter number");
    return;
  }

  Double_t plist[3] = {Double_t(_iPar), 0., 0.};
  fFitterImp.ExecuteCommand("SET LIMIT", plist, 3);
}

void
TFractionFitterSimul::FCN(Int_t& _npar, Double_t* _gin, Double_t& _f, Double_t* _par, Int_t _flag)
{
  _f = 0.;

  TObjArray& fitters(TFractionFitterSimul::fgFitters);

  TDirectory* pwd(_flag == 3 ? gDirectory : 0);

  for(Int_t iVar(0); iVar != fitters.GetEntries(); ++iVar){
    TFractionFitter* fitter(static_cast<TFractionFitter*>(fitters.At(iVar)));
    fractionFitter->SetObjectFit(fitter);

    Double_t result(0.);
    if(_flag == 3){
      TString dirName(TString::Format("TFractionFitterSimul_var%d", iVar));
      if(!pwd->cd(dirName)){
        std::cerr << "TFractionFitterSimul::FCN: Cannot cd to directory " + dirName << std::endl;
        return;
      }
    }

    TFractionFitFCN(_npar, _gin, result, _par, _flag);

    if(_flag == 3) pwd->cd();

    _f += result;
  }
}
