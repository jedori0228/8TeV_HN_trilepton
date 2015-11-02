////// eemu channel //////

#include "/Users/jskim/Documents/CMS/header/Muon.h"
#include "/Users/jskim/Documents/CMS/header/Electron.h"
#include "/Users/jskim/Documents/CMS/header/Event.h"
#include "/Users/jskim/Documents/CMS/header/trilepton.h"
#include "/Users/jskim/Documents/CMS/header/status.h"
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TStopwatch.h>
#include <TLegend.h>
#include <THStack.h>

#define W_width 2.1
#define MET_width 7
#define N_chain 24

using namespace std;

struct Fit{
  double chi2;
  double MET_fit;
  double lambda_fit;
};
double chisquare(double M_fit, double lambda);
void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, TString pm);

/////////////////////
/// main function ///
/////////////////////

void chi2_lambda_stack(){
  TStopwatch totaltime;
  totaltime.Start();
  int sample, isPU;
  
  TChain *chain_mu[N_chain];
  TChain *chain_el[N_chain];
  TChain *chain_event[N_chain];
  TChain *chain_gen[N_chain];
  
  cout << "==========================" << endl
       << "40, 50, 60 : signal" << endl
       << "==========================" << endl
       << "Run : ";
  cin >> sample;
  cout << endl;
  cout << "==========================" << endl
  << "PileUP? (0 = no, 1 = yes) : ";
  cin >> isPU;
  cout << "==========================" << endl;
  
  /// [0] = Signal ///
  chain_mu[0] = new TChain("Muon");
  chain_el[0] = new TChain("Electron");
  chain_event[0] = new TChain("Event");
  chain_mu[0]->Add("./nocut/signal/trilepton_eemu_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  chain_el[0]->Add("./nocut/signal/trilepton_eemu_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  chain_event[0]->Add("./nocut/signal/trilepton_eemu_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  
  /// [1] = Data ///
  chain_mu[1] = new TChain("Muon");
  chain_el[1] = new TChain("Electron");
  chain_event[1] = new TChain("Event");
  chain_mu[1]->Add("./nocut/data/*.root");
  chain_el[1]->Add("./nocut/data/*.root");
  chain_event[1]->Add("./nocut/data/*.root");
  
  /// [2]~ = Bkg ///
  TString chainlist[N_chain] = {"signal", "data",
                                "DY10to50", "DY50plus", "Zbb",  // DY
                                "WZtollln_mg", "WZtollqq_mg", "WZtoqqln_mg", "ZZtollll_mg", "ZZtollnn_mg", "ZZtollqq_mg", "WW_mg",  // VV
                                "Wbb", // Wjets
                                "topDIL", // Top
                                "TTG", "TTWW", "WWG", "WWW", "WWZ", "WZZ", "ZZZ", "ttZ", "HtoWW", "Wtollln"}; // others
  for(int i=2; i<N_chain; i++){
    chain_mu[i] = new TChain("Muon");
    chain_el[i] = new TChain("Electron");
    chain_event[i] = new TChain("Event");
    chain_mu[i]->Add("./nocut/bkg/trilepton_eemu_SK"+chainlist[i]+"_5_3_14.root");
    chain_el[i]->Add("./nocut/bkg/trilepton_eemu_SK"+chainlist[i]+"_5_3_14.root");
    chain_event[i]->Add("./nocut/bkg/trilepton_eemu_SK"+chainlist[i]+"_5_3_14.root");
  }
  
  Muon *mu[N_chain];
  Electron *el[N_chain];
  Event *evt[N_chain];

  double MET, METphi, nu_Pz, x, W_mass_fit, MET_fit, weight;
  int N_entry;
  struct Fit fitresult_p, fitresult_m;
  TString pm;
  Color_t fillcolors[N_chain] = {0, kRed, // signal and data. data random color
                                 kAzure-5, kAzure-5, kAzure-5,
                                 kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5,
                                 kYellow-5,
                                 kRed-5,
                                 kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kOrange-4};
  
  TH1D *hist_HN[N_chain];
  TH1D *hist_W_real[N_chain];
  TH1D *hist_l1l2[N_chain];
  TH1D *hist_l2l3[N_chain];
  
  THStack *bkgstack_HN = new THStack("bkgstack_HN", "bkgstack_HN");
  THStack *bkgstack_W_real = new THStack("bkgstack_W_real", "bkgstack_W_real");
  THStack *bkgstack_l1l2 = new THStack("bkgstack_l1l2", "bkgstack_l1l2");
  THStack *bkgstack_l2l3 = new THStack("bkgstack_l2l3", "bkgstack_l2l3");
  
  TLorentzVector lep[3], W_real, HN, nu, gen_nu;

  for(int k=0; k<N_chain; k++){
    N_entry = 0;
    mu[k] = new Muon(chain_mu[k], 1);
    el[k] = new Electron(chain_el[k], 2);
    evt[k] = new Event(chain_event[k]);
 
    hist_HN[k] = new TH1D("hist_HN_"+chainlist[k], "HN invariant mass", 150./10., 0., 150.);
    hist_W_real[k] = new TH1D("hist_W_real_"+chainlist[k], "W_real invariant mass", 30./1., 70., 100.);
    hist_l1l2[k] = new TH1D("hist_l1l2_"+chainlist[k], "l1l2 invariant mass", 200./5., 0., 200.);
    hist_l2l3[k] = new TH1D("hist_l2l3_"+chainlist[k], "l2l3 invariant mass", 200./5., 0., 200.);
    
    cout << endl << "Running Sample : " << chainlist[k] << endl;
    
  
    for(int i=0; i<chain_mu[k]->GetEntries(); i++){
      
      printeventnumber(i, chain_mu[k]->GetEntries());
      
      chain_mu[k]->GetEntry(i);
      chain_el[k]->GetEntry(i);
      chain_event[k]->GetEntry(i);
      
      if( el[k]->isMySel(0) && el[k]->isMySel(1) && mu[k]->isHNLooseMuon(0)){
      //if( mu[k]->PassID("TightMuon", 0) && mu[k]->PassID("TightMuon", 1) && mu[k]->PassID("TightMuon", 2) ){
      //if( mu[k]->isHNTightMuon(0) && mu[k]->isHNTightMuon(1) && mu[k]->isHNTightMuon(2) ){
      //if( 1 ){ // NoHNLoose cut
        lep[0].SetPxPyPzE(el[k]->Px[0],el[k]->Py[0],el[k]->Pz[0],el[k]->E[0]);
        lep[1].SetPxPyPzE(mu[k]->Px[0],mu[k]->Py[0],mu[k]->Pz[0],mu[k]->E[0]);
        lep[2].SetPxPyPzE(el[k]->Px[1],el[k]->Py[1],el[k]->Pz[1],el[k]->E[1]);

        MET = evt[k]->MET;
        METphi = evt[k]->METphi;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        
        fit(&fitresult_p, 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "p");
        fit(&fitresult_m, 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "m");
        
        pm = "p";
        /*
        if(fitresult_p.chi2 > fitresult_m.chi2){
          fitresult_p.chi2 = fitresult_m.chi2;
          fitresult_p.MET_fit = fitresult_m.MET_fit;
          fitresult_p.lambda_fit = fitresult_m.lambda_fit;
          pm = "m";
        }
        */
        
        if(fitresult_p.chi2 < 1.0 && fitresult_p.lambda_fit == 0){  // chi2 < 1 && lambda_fit == 0 cut applied
        //if(fitresult_p.chi2 < 1.0){  // chi2 < 1.0 cut applied
        //if(1){ // Loose Muon cut only
          N_entry++;
          PutNuPz(&nu, solveqdeq(80.4, lep[0]+lep[1]+lep[2], fitresult_p.MET_fit, METphi, pm));
          HN = lep[1] + lep[2] + nu;
          W_real = HN + lep[0];
          
          if(!isPU) weight = evt[k]->weight;
          if(isPU)  weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
          
          hist_HN[k]->Fill(HN.M(), weight);
          hist_W_real[k]->Fill(W_real.M(), weight);
          //hist_l1l2[k]->Fill( ( lep[0] + lep[1] ).M() , weight);
          //hist_l2l3[k]->Fill( ( lep[1] + lep[2] ).M() , weight);
          
        } // chi2 < 1 cut
      } // Loosemuon cut
    } // event loop
    
    cout << "weight = " << evt[k]->weight << endl
         << "entry = " << N_entry << endl;
    
    /*
    if(k==4){ // Add two DY process
      hist_HN[4]->Add(hist_HN[3]);
      hist_W_real[4]->Add(hist_W_real[3]);
    }
    */
    
    hist_HN[k]->SetFillColor(fillcolors[k]);
    hist_W_real[k]->SetFillColor(fillcolors[k]);
    hist_l1l2[k]->SetFillColor(fillcolors[k]);
    hist_l2l3[k]->SetFillColor(fillcolors[k]);
    
    if(k>=2){ // 24th is w_llln
      bkgstack_HN->Add(hist_HN[k]);
      bkgstack_W_real->Add(hist_W_real[k]);
      bkgstack_l1l2->Add(hist_l1l2[k]);
      bkgstack_l2l3->Add(hist_l2l3[k]);
    }
  } // chain loop
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  bkgstack_HN->Draw("hist");
  bkgstack_HN->SetMaximum(100);
  bkgstack_HN->SetTitle("HN");
  bkgstack_HN->GetXaxis()->SetTitle("m(ll#nu) [GeV]");
  bkgstack_HN->GetYaxis()->SetTitle("Event / 10 GeV");
  //hist_HN[0]->SetLineColor(kRed);
  //hist_HN[0]->Draw("histsame");
  hist_HN[1]->SetMarkerStyle(3);
  hist_HN[1]->SetMarkerSize(1);
  hist_HN[1]->Draw("APsameE1");
  hist_HN[0]->SetLineColor(kRed);
  if(sample == 40)  hist_HN[0]->Scale(0.1);
  if(sample == 60)  hist_HN[0]->Scale(0.1);
  hist_HN[0]->Draw("histsame");
  TLegend *lg1 = new TLegend(0.9, 0.6, 1, 0.9);
  lg1->AddEntry(hist_HN[1], "data", "p");
  if(sample == 40)  lg1->AddEntry(hist_HN[0], "signal*0.1", "l");
  if(sample == 60)  lg1->AddEntry(hist_HN[0], "signal*0.1", "l");
  lg1->AddEntry(hist_HN[23], "Wtollln", "f");
  lg1->AddEntry(hist_HN[14], "other", "f");
  lg1->AddEntry(hist_HN[13], "top", "f");
  lg1->AddEntry(hist_HN[12], "Wjets", "f");
  lg1->AddEntry(hist_HN[5], "VV", "f");
  lg1->AddEntry(hist_HN[2], "DY", "f");
  lg1->Draw();


  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  bkgstack_W_real->Draw("hist");
  bkgstack_W_real->SetMaximum(100);
  bkgstack_W_real->SetTitle("W");
  bkgstack_W_real->GetXaxis()->SetTitle("m(lll#nu) [GeV]");
  bkgstack_W_real->GetYaxis()->SetTitle("Event / 1 GeV");
  //hist_W_real[0]->SetLineColor(kRed);
  //hist_W_real[0]->Draw("histsame");
  hist_W_real[1]->SetMarkerStyle(3);
  hist_W_real[1]->SetMarkerSize(1);
  hist_W_real[1]->Draw("APsameE1");
  hist_W_real[0]->SetLineColor(kRed);
  if(sample == 40)  hist_W_real[0]->Scale(0.1);
  if(sample == 60)  hist_W_real[0]->Scale(0.1);
  hist_W_real[0]->Draw("histsame");
  TLegend *lg2 = new TLegend(0.9, 0.6, 1, 0.9);
  lg2->AddEntry(hist_W_real[1], "data", "p");
  if(sample == 40)  lg2->AddEntry(hist_W_real[0], "signal*0.1", "l");
  if(sample == 60)  lg2->AddEntry(hist_W_real[0], "signal*0.1", "l");
  lg2->AddEntry(hist_W_real[23], "Wtollln", "f");
  lg2->AddEntry(hist_W_real[14], "other", "f");
  lg2->AddEntry(hist_W_real[13], "top", "f");
  lg2->AddEntry(hist_W_real[12], "Wjets", "f");
  lg2->AddEntry(hist_W_real[5], "VV", "f");
  lg2->AddEntry(hist_W_real[2], "DY", "f");
  lg2->Draw();

  /*
  TCanvas *c3 = new TCanvas("c3", "", 800, 600);
  bkgstack_l1l2->Add(hist_l1l2[0]);
  bkgstack_l1l2->Draw("hist");
  bkgstack_l1l2->SetMaximum(60);
  bkgstack_l1l2->SetTitle("l_{1}l_{2} Invariant Mass");
  bkgstack_l1l2->GetXaxis()->SetTitle("m(l_{1}l_{2}) [GeV]");
  bkgstack_l1l2->GetYaxis()->SetTitle("Event / 5GeV");
  //hist_l1l2[0]->SetLineColor(kRed);
  //hist_l1l2[0]->Draw("histsame");
  hist_l1l2[1]->SetMarkerStyle(3);
  hist_l1l2[1]->SetMarkerSize(1);
  hist_l1l2[1]->Draw("APsameE1");
  TLegend *lg3 = new TLegend(0.9, 0.6, 1, 0.9);
  lg3->AddEntry(hist_l1l2[1], "data", "p");
  lg3->AddEntry(hist_l1l2[0], "signal", "f");
  lg3->AddEntry(hist_l1l2[14], "other", "f");
  lg3->AddEntry(hist_l1l2[13], "top", "f");
  lg3->AddEntry(hist_l1l2[12], "Wjets", "f");
  lg3->AddEntry(hist_l1l2[5], "VV", "f");
  lg3->AddEntry(hist_l1l2[2], "DY", "f");
  lg3->Draw();
  
  TCanvas *c4 = new TCanvas("c4", "", 800, 600);
  bkgstack_l2l3->Add(hist_l2l3[0]);
  bkgstack_l2l3->Draw("hist");
  bkgstack_l2l3->SetMaximum(60);
  bkgstack_l2l3->SetTitle("l_{2}l_{3} Invariant Mass");
  bkgstack_l2l3->GetXaxis()->SetTitle("m(l_{2}l_{3}) [GeV]");
  bkgstack_l2l3->GetYaxis()->SetTitle("Event / 5GeV");
  //hist_l2l3[0]->SetLineColor(kRed);
  //hist_l2l3[0]->Draw("histsame");
  hist_l2l3[1]->SetMarkerStyle(3);
  hist_l2l3[1]->SetMarkerSize(1);
  hist_l2l3[1]->Draw("APsameE1");
  TLegend *lg4 = new TLegend(0.9, 0.6, 1, 0.9);
  lg4->AddEntry(hist_l2l3[1], "data", "p");
  lg4->AddEntry(hist_l2l3[0], "signal", "f");
  lg4->AddEntry(hist_l2l3[14], "other", "f");
  lg4->AddEntry(hist_l2l3[13], "top", "f");
  lg4->AddEntry(hist_l2l3[12], "Wjets", "f");
  lg4->AddEntry(hist_l2l3[5], "VV", "f");
  lg4->AddEntry(hist_l2l3[2], "DY", "f");
  lg4->Draw();
  */
}


double chisquare(double M_fit, double lambda){
  return TMath::Power((M_fit-80.4)/W_width, 2) + lambda*lambda;
}

void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, TString pm){
  double dlambda = 2*range_lambda/n_lambda, lambda, x = 999999, x_temp, MET_temp; // step size
  TLorentzVector nu, W_real;

  
  for(int i=0; i<n_lambda; i++){
    lambda = -range_lambda+dlambda*i;
    MET_temp = MET_recon + lambda*7.0;
    nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
    
    if(MET_temp >= 0){ // MET should be positive
      PutNuPz(&nu, solveqdeq(80.4, l1l2l3, MET_temp, METphi_recon, pm));
      W_real = l1l2l3 + nu;
      x_temp = chisquare(W_real.M(), lambda);
    
      if(x_temp < x){
        //cout << lambda << '\t' << x_temp-lambda*lambda << endl;
        x = x_temp;
        fitresult->chi2 = x;
        fitresult->MET_fit = MET_temp;
        fitresult->lambda_fit = -range_lambda+dlambda*i;
      }
    }
  }
  
}
