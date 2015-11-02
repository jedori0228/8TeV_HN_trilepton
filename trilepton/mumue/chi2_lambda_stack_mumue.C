////// mumue channel //////

#include "/Users/jskim/Documents/CMS/header/Muon.h"
#include "/Users/jskim/Documents/CMS/header/Electron.h"
#include "/Users/jskim/Documents/CMS/header/Jet.h"
#include "/Users/jskim/Documents/CMS/header/Event.h"
#include "/Users/jskim/Documents/CMS/header/Trilepton.h"
#include "/Users/jskim/Documents/CMS/header/Status.h"
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TStopwatch.h>
#include <TLegend.h>
#include <THStack.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#define W_width 2.1
#define MET_width 7
#define N_chain 24

using namespace std;

struct Fit{
  double chi2;
  double MET_fit;
  double lambda_fit;
  double Pz_fit;
  bool isSamePz;
  Fit():  chi2(0), MET_fit(0), lambda_fit(0), Pz_fit(0), isSamePz(false){}
};
double chisquare(double M_fit, double lambda);
void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, TString pm);

/////////////////////
/// main function ///
/////////////////////

void chi2_lambda_stack_mumue(){
  TStopwatch totaltime;
  totaltime.Start();
  int sample, isPU, ischi2;
  
  TChain *chain_mu[N_chain];
  TChain *chain_el[N_chain];
  TChain *chain_jet[N_chain];
  TChain *chain_event[N_chain];
  
  cout << "==========================" << endl
  << "40, 50, 60, 100(no gen lvl info) : signal" << endl
  << "==========================" << endl
  << "Run : ";
  cin >> sample;
  cout << endl;
  
  cout << "==========================" << endl
  << "PileUP? (0 = no, 1 = yes) : ";
  cin >> isPU;
  cout << "==========================" << endl
  << "chi2 fit? (0 = no, 1 = yes) :";
  cin >> ischi2;
  cout << endl;

  /*
  /// [0] = Signal ///
  chain_mu[0] = new TChain("Muon");
  chain_el[0] = new TChain("Electron");
  chain_event[0] = new TChain("Event");
  chain_mu[0]->Add("./files/signal/trilepton_mumue_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  chain_el[0]->Add("./files/signal/trilepton_mumue_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  chain_event[0]->Add("./files/signal/trilepton_mumue_SKHN"+TString::Itoa(sample, 10)+"_lep_5_3_14.root");
  */
  
  /// [1] = Data ///
  chain_mu[1] = new TChain("Muon");
  chain_el[1] = new TChain("Electron");
  chain_event[1] = new TChain("Event");
  chain_jet[1] = new TChain("Jet");
  chain_mu[1]->Add("./files/data/*.root");
  chain_el[1]->Add("./files/data/*.root");
  chain_event[1]->Add("./files/data/*.root");
  chain_jet[1]->Add("./files/data/*.root");
  
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
    chain_jet[i] = new TChain("Jet");
    chain_mu[i]->Add("./files/bkg/trilepton_mumue_SK"+chainlist[i]+"_5_3_14.root");
    chain_el[i]->Add("./files/bkg/trilepton_mumue_SK"+chainlist[i]+"_5_3_14.root");
    chain_event[i]->Add("./files/bkg/trilepton_mumue_SK"+chainlist[i]+"_5_3_14.root");
    chain_jet[i]->Add("./files/bkg/trilepton_mumue_SK"+chainlist[i]+"_5_3_14.root");
  }
  
  Muon *mu[N_chain];
  Electron *el[N_chain];
  Event *evt[N_chain];
  Jet *jet[N_chain];
           
  double MET, METphi, nu_Pz, W_mass_fit, MET_fit, weight, PUweight_sum;
  int N_entry, HN_x_min=0, HN_x_max=100, HN_dx=10, W_x_min=70, W_x_max=120, W_dx=5, solution_selection, smaller;
  struct Fit fitresult[2];
  Color_t fillcolors[N_chain] = {0, kRed-5, // signal and data. data random color
    kAzure-5, kAzure-5, kAzure-5,
    kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5,
    kYellow-5,
    kRed-5,
    kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kOrange-4};
  double gamma_star_bin[16];
  TString gamma_star_bin_label[15];
  for(int i=0; i<15; i++){
    gamma_star_bin[i] = 0.5*i;
    if(i%2==0) gamma_star_bin_label[i] = "";
    else      gamma_star_bin_label[i] = "            "+TString::Itoa(i/2+1, 10)+".0";
    //cout << gamma_star_bin[i] << '\t' << gamma_star_bin_label[i] << endl;
  }
  gamma_star_bin[15] = 8;
  gamma_star_bin_label[14] = "        > 7.0 GeV";
  
  TH1D *hist_dimuon[N_chain];
  
  TH1D *hist_dimuon_err = new TH1D("hist_dimuon_err", "", (50.-0.)/1., 0., 50.);
  
  THStack *bkgstack_dimuon = new THStack("bkgstack_dimuon", "");
  
  TLorentzVector mu_4vec[2], el_4vec, W_real, HN, nu, gen_nu, gamma_star, z_candidate, selection_nu[2];

  for(int k=1; k<N_chain; k++){
    N_entry = 0;
    PUweight_sum=0;
    
    mu[k] = new Muon(chain_mu[k], 2);
    el[k] = new Electron(chain_el[k], 1);
    jet[k] = new Jet(chain_jet[k], 4);
    evt[k] = new Event(chain_event[k]);
 
    hist_dimuon[k] = new TH1D("hist_dimuon"+TString::Itoa(k,10), "", (50.-0.)/1., 0., 50.);

    
    cout << endl << "Running Sample : " << chainlist[k] << endl;
    
  
    for(int i=0; i<chain_mu[k]->GetEntries(); i++){
      
      printeventnumber(i, chain_mu[k]->GetEntries());
      
      chain_mu[k]->GetEntry(i);
      chain_el[k]->GetEntry(i);
      chain_jet[k]->GetEntry(i);
      chain_event[k]->GetEntry(i);
      
      if( mu[k]->isHNLooseMuon(0) && mu[k]->isHNLooseMuon(1)){
      //if( mu[k]->PassID("TightMuon", 0) && mu[k]->PassID("TightMuon", 1) && mu[k]->PassID("TightMuon", 2) ){
      //if( mu[k]->isHNTightMuon(0) && mu[k]->isHNTightMuon(1) && mu[k]->isHNTightMuon(2) ){
      //if( 1 ){ // NoHNLoose cut
        mu_4vec[0].SetPxPyPzE(mu[k]->Px[0],mu[k]->Py[0],mu[k]->Pz[0],mu[k]->E[0]);
        mu_4vec[1].SetPxPyPzE(mu[k]->Px[1],mu[k]->Py[1],mu[k]->Pz[1],mu[k]->E[1]);
        el_4vec.SetPxPyPzE(el[k]->Px[0],el[k]->Py[0],el[k]->Pz[0],el[k]->E[0]);

        MET = evt[k]->MET;
        METphi = evt[k]->METphi;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        
        if(ischi2){
          fit(&fitresult[0], 3, 1000, mu_4vec[0]+mu_4vec[1]+el_4vec[0], MET, METphi, "m");
          fit(&fitresult[1], 3, 1000, mu_4vec[0]+mu_4vec[1]+el_4vec[0], MET, METphi, "p"); // 0 = minus, 1 = plus
          PutNuPz(&selection_nu[0], fitresult[0].Pz_fit);
          PutNuPz(&selection_nu[1], fitresult[1].Pz_fit);
        }
        else{
          PutNuPz(&selection_nu[0], solveqdeq(80.4, mu_4vec[0]+mu_4vec[1]+el_4vec[0], MET, METphi, "m"));
          PutNuPz(&selection_nu[1], solveqdeq(80.4, mu_4vec[0]+mu_4vec[1]+el_4vec[0], MET, METphi, "p")); // 0 = minus, 1 = plus
        }
        
        if( selection_nu[0].Pz() == selection_nu[1].Pz() ){
          solution_selection = 0; // 0, 1 상관 없으므로
        }
        else{
          if(fabs(selection_nu[0].Pz()) > fabs(selection_nu[1].Pz())){
            smaller = 1;
          }
          else{
            smaller = 0;
          }
          
          solution_selection = smaller;
        }
        
        if(!isPU) weight = evt[k]->weight;
        else      weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
        
        if(mu[k]->Charge[0] != mu[k]->Charge[1] && evt[k]->n_bjet>=1 && jet[k]->n_Jet>=1){
        //if(mu[k]->Charge[0] != mu[k]->Charge[1] && mu[0].DeltaR(mu[1]) < 1.0){
          N_entry++;
          PUweight_sum+=evt[k]->PU_reweight;
          hist_dimuon[k]->Fill( (mu_4vec[0]+mu_4vec[1]).M(), weight );
          
        } // additional cuts
      } // Object Selection
    } // event loop
    
    cout << "weight = " << evt[k]->weight << endl
    << "entry = " << N_entry << endl
    << "PUweight = " << PUweight_sum/(double)N_entry << endl;
    
    /// apply weight after filling histogram ///
    //hist_HN[k]->Scale(weight);
    //hist_W_real[k]->Scale(weight);
    //hist_dR[k]->Scale(weight);
    //hist_gamma_star[k]->Scale(weight);
    
    /*
     if(k==4){ // Add two DY process
     hist_HN[4]->Add(hist_HN[3]);
     hist_W_real[4]->Add(hist_W_real[3]);
     }
     */
    
    hist_dimuon[k]->SetFillColor(fillcolors[k]);
    
    /// merging same-type histograms
    if(k>=15 && k<=22){
      hist_dimuon[14]->Add(hist_dimuon[k]);
    }
    if(k>=6 && k<=11){
      hist_dimuon[5]->Add(hist_dimuon[k]);
    }
    if(k==3 || k==4){
      hist_dimuon[2]->Add(hist_dimuon[k]);
    }
    if(k==14 || k==23 || k==5 || k==12 || k==13 || k==2){
      bkgstack_dimuon->Add(hist_dimuon[k]);
    }
    if(k>=2){
      hist_dimuon_err->Add(hist_dimuon[k]);
    }
  } // chain loop
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  /// dimuon ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->cd();
  ///bkg///
  bkgstack_dimuon->Draw("hist");
  bkgstack_dimuon->SetMaximum(100); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_dimuon->SetTitle("Dimuon(OS) Invariant Mass");
  //bkgstack_dimuon->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  bkgstack_dimuon->GetYaxis()->SetTitle("Event / 1 GeV");
  ///data///
  hist_dimuon[1]->SetMarkerStyle(3);
  hist_dimuon[1]->SetMarkerSize(1);
  hist_dimuon[1]->Draw("APsameE1");
  ///signal///
  //hist_dimuon[0]->SetLineColor(kRed);
  //hist_dimuon[0]->Draw("histsame");
  ///err///
  //hist_dimuon_err->SetFillStyle(3004);
  //hist_dimuon_err->SetFillColor(kBlue);
  //hist_dimuon_err->Draw("sameE2");
  ///legend///
  TLegend *lg1 = new TLegend(0.65, 0.6, 0.9, 0.9);
  lg1->AddEntry(hist_dimuon[1], "data", "p");
  //lg1->AddEntry(hist_dimuon[0], "HN"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  lg1->AddEntry(hist_dimuon[2], "DY", "f");
  lg1->AddEntry(hist_dimuon[13], "top", "f");
  lg1->AddEntry(hist_dimuon[5], "VV", "f");
  lg1->AddEntry(hist_dimuon[23], "Wtollln", "f");
  lg1->AddEntry(hist_dimuon[12], "Wjets", "f");
  lg1->AddEntry(hist_dimuon[14], "other", "f");
  lg1->Draw();
}


double chisquare(double M_fit, double lambda){
  return TMath::Power((M_fit-80.4)/W_width, 2) + lambda*lambda;
}

void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, TString pm){
  double dlambda = range_lambda/(n_lambda-1), lambda, x, x_temp, MET_temp, Pz_temp; // step size
  TLorentzVector nu, W_real;
  
  lambda = 0;
  MET_temp = MET_recon;
  nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
  Pz_temp = solveqdeq(&(fitresult->isSamePz), 80.4, l1l2l3, MET_temp, METphi_recon, pm);
  PutNuPz(&nu, Pz_temp);
  W_real = l1l2l3 + nu;
  x_temp = chisquare(W_real.M(), lambda);
  x = x_temp;
  fitresult->chi2 = x;
  fitresult->MET_fit = MET_temp;
  fitresult->lambda_fit = lambda;
  fitresult->Pz_fit = Pz_temp;
  //cout << x_temp << '\t' << W_real.M() << '\t' << fitresult->isSamePz << endl;
  
  for(int i=0; i<2*(n_lambda-1); i++){
    if(lambda==0){
      lambda = dlambda;
    }
    else if(lambda > 0){
      lambda = -lambda;
    }
    else{ // lambda < 0
      lambda = -lambda+dlambda;
    }
    //cout << lambda << endl;
    
    MET_temp = MET_recon + lambda*7.0;
    nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
    
    if(MET_temp >= 0){ // MET should be positive
      Pz_temp = solveqdeq(&(fitresult->isSamePz), 80.4, l1l2l3, MET_temp, METphi_recon, pm);
      PutNuPz(&nu, Pz_temp);
      W_real = l1l2l3 + nu;
      x_temp = chisquare(W_real.M(), lambda);
      if(W_real.M()==80.4) break; // lambda = 0 solvable
      if(x_temp < x){
        //cout << "old chi2 : " << x << '\t' << "new chi2 : " << x_temp << endl;
        x = x_temp;
        fitresult->chi2 = x;
        fitresult->MET_fit = MET_temp;
        fitresult->lambda_fit = -range_lambda+dlambda*i;
        fitresult->Pz_fit = Pz_temp;
      }
    }
  }
  
}

