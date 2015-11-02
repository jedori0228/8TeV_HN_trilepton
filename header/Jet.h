#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <TMath.h>

#define MaxJetNumber 10

using namespace std;

class Jet{
public:
  double Px[MaxJetNumber];
  double Py[MaxJetNumber];
  double Pz[MaxJetNumber];
  double E[MaxJetNumber];
  double Eta[MaxJetNumber];
  double Pt[MaxJetNumber];
  double PFJetTrackCountingHighPurBTag[MaxJetNumber];
  double BtagProb[MaxJetNumber];
  double CombinedSecVertexBtag[MaxJetNumber];
  double n_Jet;
  
  Jet(TChain *chain, int n_jetini);
};

Jet::Jet(TChain *chain, int n_jetini){
  TString varname[9] = {"Px", "Py", "Pz", "E", "Eta", "Pt", "PFJetTrackCountingHighPurBTag", "BtagProb", "CombinedSecVertexBtag"};
  chain->SetBranchAddress("n_Jet", &n_Jet);
  for(int i=0; i<n_jetini; i++){
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[0], Px+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[1], Py+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[2], Pz+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[3], E+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[4], Eta+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[5], Pt+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[6], PFJetTrackCountingHighPurBTag+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[7], BtagProb+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[8], CombinedSecVertexBtag+i);
  }
}





