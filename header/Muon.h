#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <TMath.h>

#define MaxMuonNumber 3

using namespace std;

class Muon{
public:
  double Px[MaxMuonNumber];
  double Py[MaxMuonNumber];
  double Pz[MaxMuonNumber];
  double E[MaxMuonNumber];
  double LeptonRelIso[MaxMuonNumber];
  double Eta[MaxMuonNumber];
  double dXY[MaxMuonNumber];
  double dZ[MaxMuonNumber];
  double GlobalChi2[MaxMuonNumber];
  double IsoHcalVeto[MaxMuonNumber];
  double IsoEcalVeto[MaxMuonNumber];
  double IsGlobal[MaxMuonNumber];
  double validHits[MaxMuonNumber];
  double validPixHits[MaxMuonNumber];
  double validStations[MaxMuonNumber];
  double ActiveLayer[MaxMuonNumber];
  double IsPf[MaxMuonNumber];
  double IsTracker[MaxMuonNumber];
  double Charge[MaxMuonNumber];
  double Pt[MaxMuonNumber];
  double n_LooseMuon;
  
  Muon(TChain *chain, int n_muonini);
  bool PassID(TString ID, int i);
  bool isLooseMuon(int i);
  bool isHNLooseMuon(int i, double min_pt);
  bool isHNTightMuon(int i, double min_pt);
};

Muon::Muon(TChain *chain, int n_muonini){
  TString varname[20] = {"Px", "Py", "Pz", "E", "LeptonRelIso", "Eta", "dXY", "dZ", "GlobalChi2", "IsoHcalVeto", "IsoEcalVeto", "IsGlobal", "validHits", "validPixHits", "validStations", "ActiveLayer", "IsPf", "IsTracker", "Charge", "Pt"};
  chain->SetBranchAddress("n_LooseMuon", &n_LooseMuon);
  for(int i=0; i<n_muonini; i++){
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[0], Px+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[1], Py+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[2], Pz+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[3], E+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[4], LeptonRelIso+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[5], Eta+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[6], dXY+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[7], dZ+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[8], GlobalChi2+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[9], IsoHcalVeto+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[10], IsoEcalVeto+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[11], IsGlobal+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[12], validHits+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[13], validPixHits+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[14], validStations+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[15], ActiveLayer+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[16], IsPf+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[17], IsTracker+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[18], Charge+i);
    chain->SetBranchAddress(TString::Itoa(i, 10)+"_"+varname[19], Pt+i);
  }
}

bool Muon::PassID(TString id, int i){
  if(id=="LooseMuon"){
    if(  !(IsPf[i] == 1)  ){
      return false;
    }
    if(  !(IsGlobal[i] == 1 || IsTracker[i] == 1)  ){
      return false;
    }
    return true;
  }
  
  else if(id=="TightMuon"){
    if(  !(IsPf[i] == 1)  ){
      return false;
    }
    if(  !(IsGlobal[i] == 1)  ){
      return false;
    }
    if(  validHits[i] == 0  ){
      return false;
    }
    if(  validPixHits[i] == 0  ){
      return false;
    }
    if(  validStations[i] <= 1  ){
      return false;
    }
    if(  ActiveLayer[i] <= 5  ){
      return false;
    }
    if(  fabs(dXY[i]) >= 0.2  ){
      return false;
    }
    if(  fabs(dZ[i]) >= 0.5  ){
      return false;
    }
    if(  GlobalChi2[i] >= 10.  ){
      return false;
    }
    return true;
  }
  
  else return true;
}

bool Muon::isHNLooseMuon(int i, double min_pt){
  if(  Pt[i] < min_pt  )                             return false;
  if(  !(fabs(Eta[i]) < 2.4)  )                      return false;
  if(  !PassID("LooseMuon", i)  )                    return false;
  if(  !(LeptonRelIso[i] < 0.1)  )                   return false;
  if(  !(IsGlobal[i] == 1)  )                        return false;
  if(  validHits[i] == 0  )                          return false;
  if(  validPixHits[i] == 0  )                       return false;
  if(  validStations[i] <= 1  )                      return false;
  if(  ActiveLayer[i] <= 5  )                        return false;
  if(  fabs(dXY[i]) >= 0.2  )                        return false;
  if(  fabs(dZ[i]) >= 0.1  )                         return false;
  if(  !(GlobalChi2[i] < 50.)  )                     return false;
  if(  !(IsoHcalVeto[i] < 6.0)  )                    return false;
  if(  !(IsoEcalVeto[i] < 4.0)  )                    return false;
  
  return true;
}

bool Muon::isHNTightMuon(int i, double min_pt){
  if(  Pt[i] < min_pt  ) return false;
  if(  !(fabs(Eta[i]) < 2.4)  )                      return false;
  if(  !(LeptonRelIso[i] < 0.05)  )                  return false;
  if(  !(fabs(dZ[i]) < 0.10)  )                      return false;
  if(  !(fabs(dXY[i]) < 0.005)  )                    return false;
  if(  !PassID("TightMuon", i)  )                    return false;
  if(  !(IsoHcalVeto[i] < 6.0)  )                    return false;
  if(  !(IsoEcalVeto[i] < 4.0)  )                    return false;
  
  return true;
}






