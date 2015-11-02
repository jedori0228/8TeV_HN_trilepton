#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TChain.h>
#include "trilepton.h"
#include "Selection.h"

void maincode(){
  
  TChain *chain = new TChain("myntp");
  chain->Add("*.root");
  

  
  TLorentzVector lep[3], nu, HN, W_real; /// recon
  TH1D *hist_HN_plus = new TH1D("hist_HN_plus", "#mue#nu_plus", 200./5. ,0., 200.);
  TH1D *hist_HN_minus = new TH1D("hist_HN_minus", "#mue#nu_plus", 200./5. ,0., 200.);
  TH1D *hist_HN_inverse_plus = new TH1D("hist_HN_inverse_plus", "#mue#nu_inverse_plus", 200./5. ,0., 200.);
  TH1D *hist_HN_inverse_minus = new TH1D("hist_HN_inverse_minus", "#mue#nu_inverse_plus", 200./5. ,0., 200.);
  TH1D *hist_W_real_plus = new TH1D("hist_W_real_plus", "e#mue#nu_plus", 200./5. ,0., 200.);
  TH1D *hist_W_real_minus = new TH1D("hist_W_real_minus", "e#mue#nu_minus", 200./5. ,0., 200.);

  
  
  for(int i=0;i<recon.GetEntries();i++){ // i-th event starts
    recon.GetEntry(i);
    
    /// TLorentzVector settings ///
    TLorentzInitialize(lep, &nu, recon);
    
    /// Solve Quadratic Euqation ///
    bool isNegative;
    double PzSol[2];
    solveqdeq(&isNegative, PzSol, lep, recon.var[12], recon.var[13]);
    
    SetNeutrinoPz(&nu, PzSol[0]);
    HN = lep[1] + lep[2] + nu;
    W_real = HN + lep[0];
    hist_HN_plus->Fill(HN.M());
    hist_W_real_plus->Fill(W_real.M());
    HN = lep[1] + lep[0] + nu;
    hist_HN_inverse_plus->Fill(HN.M());

    SetNeutrinoPz(&nu, PzSol[1]);
    HN = lep[1] + lep[2] + nu;
    W_real = HN + lep[0];
    hist_HN_minus->Fill(HN.M());
    hist_W_real_minus->Fill(W_real.M());
    HN = lep[1] + lep[0] + nu;
    hist_HN_inverse_minus->Fill(HN.M());
  

  }
  

  
  TCanvas *c1 = new TCanvas("c1", "", 1600, 1200);
  c1->Divide(2,2);
  c1->cd(1);
  hist_HN_plus->Draw("hist");
  hist_HN_plus->SetXTitle("M(#mue#nu)");
  c1->cd(2);
  hist_HN_minus->Draw("hist");
  hist_HN_minus->SetXTitle("M(#mue#nu)");
  c1->cd(3);
  hist_HN_inverse_plus->Draw("hist");
  hist_HN_inverse_plus->SetXTitle("M(#mue#nu)");
  c1->cd(4);
  hist_HN_inverse_minus->Draw("hist");
  hist_HN_inverse_minus->SetXTitle("M(#mue#nu)");

  TCanvas *c2 = new TCanvas("c2", "", 1600, 600);
  c2->Divide(2,1);
  c2->cd(1);
  hist_W_real_plus->Draw("hist");
  hist_W_real_plus->SetXTitle("M(e#mue#nu)");
  c2->cd(2);
  hist_W_real_minus->Draw("hist");
  hist_W_real_minus->SetXTitle("M(e#mue#nu)");



  
}

