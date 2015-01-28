#ifndef TagProbeTree_h
#define TagProbeTree_h

#include "TString.h"
#include "TTree.h"

class RooArgSet;
class RooCategory;
class RooRealVar;

#include <map>
#include <vector>
#include <set>

class TagProbeTreeEntry {
 public:
  TagProbeTreeEntry(unsigned char, bool, bool);
  ~TagProbeTreeEntry() {}
  void push_back(std::map<TString, bool> const&, std::map<TString, float> const&);
  void setVar(TString const&, float);
  void fill(TTree&);
  void init(unsigned, unsigned);
  void book(TTree&);
  void set(TTree&);
  double getMass(bool = false, std::set<int> const* = 0);
  unsigned const& getRunNumber() const { return runNumber_; }
  unsigned const& getEventNumber() const { return eventNumber_; }
  int const& getFlag(TString const&) const;
  float const& getVar(TString const&) const;
  unsigned getNFlags() const { return flagNames_.size(); }
  unsigned getNObjectVars() const { return objectVarNames_.size(); }
  unsigned getNCombinationVars() const { return combinationVarNames_.size(); }
  RooArgSet* getListOfObservables() const;

 private:
  unsigned char const N_;
  std::vector<TString> flagNames_;
  std::vector<TString> objectVarNames_;
  std::vector<TString> combinationVarNames_;

  unsigned char iProd_;
  unsigned runNumber_;
  unsigned eventNumber_;
  std::map<TString, int> flags_;
  std::map<TString, float> vars_;
};

#endif
