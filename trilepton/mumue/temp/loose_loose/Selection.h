#include <TMath.h>
#include <iostream>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>

/////////////////
/// Cut class ///
/////////////////
class Cut{
public:
  double var[42];
  TTree *tree;
  
  /// constructor ///
  Cut(const char* filename ){
    /// variables ///
    TFile *rootfile = new TFile(filename);
    tree = (TTree*)rootfile->Get("myntp");
    TString varname[39] = {"l_1_Pt", "l_1_RelIso", "l_1_Eta", "l_1_dXY", "l_1_dZ", "l_1_GlobalChi2", "l_1_IsoHcalVeto", "l_1_IsoEcalVeto", "l_1_IsGlobal", "l_1_validHits", "l_1_validPixHits", "l_1_validStations", "l_1_ActiveLayer", "l_2_Pt", "l_2_RelIso", "l_2_Eta", "l_2_dXY", "l_2_dZ", "l_2_GlobalChi2", "l_2_IsoHcalVeto", "l_2_IsoEcalVeto", "l_2_IsGlobal", "l_2_validHits", "l_2_validPixHits", "l_2_validStations", "l_2_ActiveLayer", "l_3_Pt", "l_3_RelIso", "l_3_Eta", "l_3_dXY", "l_3_dZ", "l_3_GlobalChi2", "l_3_IsoHcalVeto", "l_3_IsoEcalVeto", "l_3_IsGlobal", "l_3_validHits", "l_3_validPixHits", "l_3_validStations", "l_3_ActiveLayer"};
    
    /// branch settings ///
    tree->SetBranchStatus("*", 0);
    for(int i=0; i<39; i++){
      tree->SetBranchStatus(varname[i], 1);
      tree->SetBranchAddress(varname[i], var+i);
    }
  }
  
  bool PtCut(double minPt);
  bool RelIsoCut(double maxRelIso);
  bool EtaCut(double maxEta);
  bool dXYCut(double maxdXY);
  bool dZCut(double maxdZ);
  bool GlobalChi2Cut(double maxGlobalChi2);
  bool IsoHcalVetoCut(double maxIsoHcal);
  bool IsoEcalVetoCut(double maxIsoEcal);
  bool IsGlobalCut();
  bool validHitsCut(int minvalidHits);
  bool validPixHitsCut(int minvalidPixHits);
  bool validStationsCut(int minvalidStations);
  bool ActiveLayerCut(int minActiveLayer);
  
  bool LooseMuon();
  bool cut_by_int(int N);
};

bool Cut::PtCut(double minPt){
  for(int i=0; i<3; i++){
    if(var[0+i*13] < minPt){
      return false;
    }
  }
  return true;
}

bool Cut::RelIsoCut(double maxRelIso){
  for(int i=0; i<3; i++){
    if(var[1+i*13] > maxRelIso){
      return false;
    }
  }
  return true;
}

bool Cut::EtaCut(double maxEta){
  for(int i=0; i<3; i++){
    if(fabs(var[2+i*13]) > maxEta){
      return false;
    }
  }
  return true;
}

bool Cut::dXYCut(double maxdXY){
  for(int i=0; i<3; i++){
    if(var[3+i*13] > maxdXY){
      return false;
    }
  }
  return true;
}

bool Cut::dZCut(double maxdZ){
  for(int i=0; i<3; i++){
    if(var[4+i*13] > maxdZ){
      return false;
    }
  }
  return true;
}

bool Cut::GlobalChi2Cut(double maxGlobalChi2){
  for(int i=0; i<3; i++){
    if(var[5+i*13] > maxGlobalChi2){
      return false;
    }
  }
  return true;
}

bool Cut::IsoHcalVetoCut(double maxIsoHcal){
  for(int i=0; i<3; i++){
    if(var[6+i*13] > maxIsoHcal){
      return false;
    }
  }
  return true;
}

bool Cut::IsoEcalVetoCut(double maxIsoEcal){
  for(int i=0; i<3; i++){
    if(var[7+i*13] > maxIsoEcal){
      return false;
    }
  }
  return true;
}

bool Cut::IsGlobalCut(){
  if(var[8+0*13]*var[8+1*13]*var[8+2*13] == 0){
    return false;
  }
  return true;
}

bool Cut::validHitsCut(int minvalidHits){
  for(int i=0; i<3; i++){
    if(var[9+i*13] < minvalidHits){
      return false;
    }
  }
  return true;
}

bool Cut::validPixHitsCut(int minvalidPixHits){
  for(int i=0; i<3; i++){
    if(var[10+i*13] < minvalidPixHits){
      return false;
    }
  }
  return true;
}

bool Cut::validStationsCut(int minvalidStations){
  for(int i=0; i<3; i++){
    if(var[11+i*13] < minvalidStations){
      return false;
    }
  }
  return true;
}

bool Cut::ActiveLayerCut(int minActiveLayer){
  for(int i=0; i<3; i++){
    if(var[12+i*13] < minActiveLayer){
      return false;
    }
  }
  return true;
}

////////////////////////////
/// Loose Muon Selection ///
////////////////////////////
bool Cut::LooseMuon(){
  bool pass;
  pass = PtCut(10.0);
  pass *= RelIsoCut(0.1);
  pass *= EtaCut(2.4);
  pass *= dXYCut(0.2);
  pass *= dZCut(0.1);
  pass *= GlobalChi2Cut(10.0);
  pass *= IsoHcalVetoCut(6.0);
  pass *= IsoEcalVetoCut(4.0);
  pass *= IsGlobalCut();
  pass *= validHitsCut(1);
  pass *= validPixHitsCut(1);
  pass *= validStationsCut(2);
  pass *= ActiveLayerCut(6);
  return pass;
}

bool Cut::cut_by_int(int N){
  bool pass;
  if(N==0)      pass = PtCut(10.0);
  else if(N==1) pass = RelIsoCut(0.1);
  else if(N==2) pass = EtaCut(2.4);
  else if(N==3) pass = dXYCut(0.2);
  else if(N==4) pass = dZCut(0.1);
  else if(N==5) pass = GlobalChi2Cut(10.0);
  else if(N==6) pass = IsoHcalVetoCut(6.0);
  else if(N==7) pass = IsoEcalVetoCut(4.0);
  else if(N==8) pass = IsGlobalCut();
  else if(N==9) pass = validHitsCut(1);
  else if(N==10) pass = validPixHitsCut(1);
  else if(N==11) pass = validStationsCut(2);
  else if(N==12) pass = ActiveLayerCut(6);
  return pass;
}


