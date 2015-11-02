#include <TChain.h>
#include <TMath.h>

#define MaxElectronNumber 3

class Electron{
public:
  double HasMatchedConvPhot[MaxElectronNumber];
  double MissingHits[MaxElectronNumber];
  double LeptonRelIsoDR03[MaxElectronNumber];
  double GsfCtfScPixChargeConsistency[MaxElectronNumber];
  double Eta[MaxElectronNumber];
  double Px[MaxElectronNumber];
  double Py[MaxElectronNumber];
  double Pz[MaxElectronNumber];
  double E[MaxElectronNumber];
  double dxy[MaxElectronNumber];
  double id_scalefactor[MaxElectronNumber];
  double Charge[MaxElectronNumber];
  double Pt[MaxElectronNumber];
  double n_LooseElectron;
  
  Electron(TChain *chain, int n_elini);
  bool isMySel(int i, double min_pt);
};

Electron::Electron(TChain *chain, int n_elini){
  TString varname[13] = {"HasMatchedConvPhot", "MissingHits", "LeptonRelIsoDR03", "GsfCtfScPixChargeConsistency", "Eta", "Px", "Py", "Pz", "E", "dxy", "id_scalefactor", "Charge", "Pt"};
  
  chain->SetBranchAddress("n_LooseElectron", &n_LooseElectron);
  for(int i=0; i<n_elini; i++){
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[0], HasMatchedConvPhot+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[1], MissingHits+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[2], LeptonRelIsoDR03+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[3], GsfCtfScPixChargeConsistency+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[4], Eta+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[5], Px+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[6], Py+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[7], Pz+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[8], E+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[9], dxy+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[10], id_scalefactor+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[11], Charge+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[12], Pt+i);
  }
}

bool Electron::isMySel(int i, double min_pt){
  if(  Pt[i] < min_pt  )                             return false;

  return true;
}