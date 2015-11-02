#include <TChain.h>

class Event{
public:
  double MET;
  double METphi;
  double n_bjet;
  double weight;
  double PU_reweight;
  
  Event(TChain *chain);
};

Event::Event(TChain *chain){
  TString varname[5] = {"MET", "METphi", "n_bjet", "weight", "PU_reweight"};
  chain->SetBranchAddress(varname[0], &MET);
  chain->SetBranchAddress(varname[1], &METphi);
  chain->SetBranchAddress(varname[2], &n_bjet);
  chain->SetBranchAddress(varname[3], &weight);
  chain->SetBranchAddress(varname[4], &PU_reweight);
  
}
