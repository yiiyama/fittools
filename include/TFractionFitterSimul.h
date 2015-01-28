#ifndef TFractionFitterSimul_H
#define TFractionFitterSimul_H

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

#ifndef ROOT_TH1
#include "TH1.h"
#endif

#ifndef ROOT_TFitter
#include "TFitter.h"
#endif

#ifndef ROOT_TFractionFitter
#include "TFractionFitter.h"
#endif

class TFractionFitterExt : public TFractionFitter {
 public:
  TFractionFitterExt();
  TFractionFitterExt(TH1*, TObjArray*, Option_t* = "");
  ~TFractionFitterExt();

  void ComputeChisquareLambda();
  Bool_t GetFitDone() const { return fFitDone; }
  TH1* GetPlot() const;
  void GetRanges(Int_t&, Int_t&, Int_t&, Int_t&, Int_t&, Int_t&) const;
  Bool_t IsExcluded(Int_t) const;
  void SetFitDone(Bool_t _done = kTRUE) { fFitDone = _done; }

 public:
  ClassDef(TFractionFitterExt, 1)
};

class TFractionFitterSimul : public TObject {
 public:
  TFractionFitterSimul();
  TFractionFitterSimul(TH1*, TObjArray*, Option_t* = "");
  TFractionFitterSimul(TObjArray const&, TObjArray const&, Option_t* = "");
  virtual ~TFractionFitterSimul();

  virtual void AddDataAndMCs(TH1*, TObjArray*);
  virtual void AddDataAndMCs(TObjArray const&, TObjArray const&);
  virtual void Constrain(Int_t, Double_t, Double_t);
  virtual void ErrorAnalysis(Double_t);
  virtual Int_t Fit(Option_t* = "");
  virtual Double_t GetChisquare() const;
  virtual Bool_t GetFitDone() const;
  virtual TFractionFitterExt* GetFractionFitter(Int_t _iVar) const { return static_cast<TFractionFitterExt*>(fFitters.At(_iVar)); }
  virtual TObjArray GetMCPredictions(Int_t) const;
  virtual Int_t GetNDF() const;
  virtual TObjArray GetPlots() const;
  virtual Double_t GetProb() const;
  virtual void GetResult(Int_t, Double_t&, Double_t&) const;
  virtual void UnConstrain(Int_t);

  static TObjArray fgFitters;
  static void FCN(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

 protected:
  Int_t fNpar;
  TObjArray fFitters;
  TFitter fFitterImp;
  TObjArray fDirectories;

 public:
  ClassDef(TFractionFitterSimul, 1)
};

#endif
