#include <TLorentzVector.h>
#include <TMath.h>
#include <iostream>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>

///////////////////
/// Recon class ///
///////////////////
class Recon_PtOrder{
public:
  TChain *chain;
  double var[14];
  
  /// constructor ///
  Recon_PtOrder(TChain *chainptr){
    chain = chainptr;
    /// variables ///
    TString varname[14] = {"mu_1_Px", "mu_1_Py", "mu_1_Pz", "mu_1_E", "el_Px", "el_Py", "el_Pz", "el_E", "mu_2_Px", "mu_2_Py", "mu_2_Pz", "mu_2_E", "MET", "METphi"};
    
    /// branch settings ///
    chain->SetBranchStatus("*", 0);
    for(int i=0;i<14;i++){
      chain->SetBranchStatus(varname[i], 1);
      chain->SetBranchAddress(varname[i], var+i);
    }
  }
  
  void GetEntry(Int_t i){
    chain->GetEntry(i);
  }
  
  Long64_t GetEntries(){
    return chain->GetEntries();
  }
  
};

/*
/////////////////
/// Gen class ///
/////////////////
class Gen{
public:
  double var[24];
  TTree *tree;
  
  /// constructor ///
  Gen(const char* filename ){
    /// variables ///
    TFile *rootfile = new TFile(filename);
    tree = (TTree*)rootfile->Get("myntp");
    TString varname[24] = {"gen_l_1_Px", "gen_l_1_Py", "gen_l_1_Pz", "gen_l_1_E", "gen_l_2_Px", "gen_l_2_Py", "gen_l_2_Pz", "gen_l_2_E", "gen_l_3_Px", "gen_l_3_Py", "gen_l_3_Pz", "gen_l_3_E", "gen_nu_Px", "gen_nu_Py", "gen_nu_Pz", "gen_nu_E", "gen_W_real_Px", "gen_W_real_Py", "gen_W_real_Pz", "gen_W_real_E", "gen_HN_Px", "gen_HN_Py", "gen_HN_Pz", "gen_HN_E"};
    
    /// branch settings ///
    tree->SetBranchStatus("*", 0);
    for(int i=0;i<24;i++){
      tree->SetBranchStatus(varname[i], 1);
      tree->SetBranchAddress(varname[i], var+i);
    }
  }
  
};
*/

////////////////////////////////////
/// Initilaizing lep, nu - recon ///
////////////////////////////////////
void TLorentzInitialize(TLorentzVector *lep, TLorentzVector *nu, Recon_PtOrder recon){
  for(int i=0;i<3;i++){
    lep[i].SetPxPyPzE(recon.var[0+4*i], recon.var[1+4*i],recon.var[2+4*i], recon.var[3+4*i]);
  }
  (*nu).SetPxPyPzE(recon.var[12]*cos(recon.var[13]), recon.var[12]*sin(recon.var[13]), 0, recon.var[12]);
}

/*
//////////////////////////////////
/// Initilaizing lep, nu - gen ///
//////////////////////////////////
void TLorentzInitialize(TLorentzVector *gen_lep, TLorentzVector *gen_nu, TLorentzVector *gen_W_real, TLorentzVector *gen_HN, Gen gen){
  for(int i=0;i<3;i++){
    gen_lep[i].SetPxPyPzE(gen.var[0+4*i], gen.var[1+4*i],gen.var[2+4*i], gen.var[3+4*i]);
  }
  (*gen_nu).SetPxPyPzE(gen.var[12], gen.var[13], gen.var[14], gen.var[15]);
  (*gen_W_real).SetPxPyPzE(gen.var[16], gen.var[17], gen.var[18], gen.var[19]);
  (*gen_HN).SetPxPyPzE(gen.var[20], gen.var[21], gen.var[22], gen.var[23]);
}
*/

////////////////////////////////
/// Solve quadratic equation ///
////////////////////////////////
void solveqdeq(bool *isNegative, double *PzSol, TLorentzVector *lep, double MET, double METphi){
  TLorentzVector l1l2l3 = lep[0] + lep[1] + lep[2];
  TLorentzVector met;
  met.SetPxPyPzE(MET*cos(METphi),
                 MET*sin(METphi),
                 0,
                 MET);

  double W_mass = 80.385;
  Double_t d = (W_mass*W_mass)-(l1l2l3.M())*(l1l2l3.M())+2.0*l1l2l3.Px()*met.Px()+2.0*l1l2l3.Py()*met.Py();
  Double_t a = l1l2l3.E()*l1l2l3.E() - l1l2l3.Pz()*l1l2l3.Pz();
  Double_t b = d*l1l2l3.Pz();
  Double_t c = l1l2l3.E()*l1l2l3.E()*met.E()*met.E()-d*d/4.0;
  if(b*b-4*a*c<0){
    *isNegative = true;
  }
  else{
    *isNegative = false;
  }
  *PzSol = (b+TMath::Sqrt(b*b-4*a*c))/(2*a);
  *(PzSol+1) = (b-TMath::Sqrt(b*b-4*a*c))/(2*a);
}

/////////////////////////////////////////////////////////
/// Solve quadratic equation with spcified W_real mass///
/////////////////////////////////////////////////////////
void solveqdeq(bool *isNegative, double *PzSol, TLorentzVector *lep, double MET, double METphi, double W_mass){
  TLorentzVector l1l2l3 = lep[0] + lep[1] + lep[2];
  TLorentzVector met;
  met.SetPxPyPzE(MET*cos(METphi),
                 MET*sin(METphi),
                 0,
                 MET);
  
  Double_t d = (W_mass*W_mass)-(l1l2l3.M())*(l1l2l3.M())+2.0*l1l2l3.Px()*met.Px()+2.0*l1l2l3.Py()*met.Py();
  Double_t a = l1l2l3.E()*l1l2l3.E() - l1l2l3.Pz()*l1l2l3.Pz();
  Double_t b = d*l1l2l3.Pz();
  Double_t c = l1l2l3.E()*l1l2l3.E()*met.E()*met.E()-d*d/4.0;
  if(b*b-4*a*c<0){
    *isNegative = true;
  }
  else{
    *isNegative = false;
  }
  *PzSol = (b+TMath::Sqrt(b*b-4*a*c))/(2*a);
  *(PzSol+1) = (b-TMath::Sqrt(b*b-4*a*c))/(2*a);
}



/////////////////////
/// SetNeutrinoPz ///
/////////////////////
void SetNeutrinoPz(TLorentzVector *nu, double Pz){
  double MET = TMath::Sqrt((*nu).Px()*(*nu).Px() + (*nu).Py()*(*nu).Py());
  (*nu).SetPxPyPzE((*nu).Px(), (*nu).Py(), Pz, TMath::Sqrt(MET*MET+Pz*Pz));

}

/////////////////////////
/// Angle in XY plane ///
/////////////////////////
double XYAngle(TLorentzVector a, TLorentzVector b){
  TLorentzVector a_temp, b_temp;
  a_temp.SetPxPyPzE(a.Px(), a.Py(), 0, 0);
  b_temp.SetPxPyPzE(b.Px(), b.Py(), 0, 0);
  return a_temp.Angle(b_temp.Vect());
}

