////// mumumu channel //////

#include "/home/jskim/Documents/CMS/header/Muon.h"
#include "/home/jskim/Documents/CMS/header/Event.h"
#include "/home/jskim/Documents/CMS/header/Jet.h"
#include "/home/jskim/Documents/CMS/header/Trilepton.h"
#include "/home/jskim/Documents/CMS/header/Status.h"
#include "/home/jskim/Documents/CMS/header/Logger.h"
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
#include <TText.h>

#define W_width 2.1
#define MET_width 7.0
#define N_chain 25
#define N_plots 5

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
int findleadlep(TLorentzVector *lep);

/////////////////////
/// main function ///
/////////////////////

void control_check_mumumu(){
  
  TH1::SetDefaultSumw2(true);
  
  TStopwatch totaltime;
  totaltime.Start();
  int sample, isPU, ischi2, additional_cut, isSave=0, isClose=0, y_max[N_plots];
  double met_x_min, met_x_max, met_dx, lead_pt_x_min, lead_pt_x_max, lead_pt_dx, lead_eta_x_min, lead_eta_x_max, lead_eta_dx, n_jets_x_min, n_jets_x_max, n_jets_dx, n_bjets_x_min, n_bjets_x_max, n_bjets_dx;
  TString savepath = "/home/jskim/Documents/CMS/trilepton/mumumu/plots/control/", muonidsel="loose";
  bool selection;
  TChain *chain_mu[N_chain], *chain_event[N_chain], *chain_jet[N_chain];

  cout
  << "===============================================" << endl
  << "40, 50, 60, 100(no gen lvl info) : signal" << endl
  << "===============================================" << endl
  << "Run : ";
  cin >> sample;
  cout << endl;
  
  cout
  << "===============================================" << endl
  << "PileUP? (0 = no, 1 = yes) : ";
  cin >> isPU;
  cout
  << "===============================================" << endl
  << "Additional Cut?" << endl
  << "  0 : No Additional Cut" << endl
  << "  1 : delta R > 0.5" << endl
  << "  2 : delta R > 0.5 & W_reco < 100 GeV" << endl
  << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
  << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda = 0" << endl
  << "  5 : delta R > 0.5 & # of b-tagged jet >= 1" << endl
  << "===============================================" << endl
  << "Cut : ";
  cin >> additional_cut;
  cout << endl;
  
  if(additional_cut==0){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=1000;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=1000;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=1000;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "before_dR/";
    ischi2 = 0;
  }
  if(additional_cut==1){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=800;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=800;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=1000;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "dR/";
    ischi2 = 0;
  }
  if(additional_cut==2){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=100;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "dR_W/";
    ischi2 = 0;
  }
  if(additional_cut==3){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "dR_chi2/";
    ischi2 = 1;
  }
  if(additional_cut==4){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "dR_chi2_lambda/";
    ischi2 = 1;
  }
  if(additional_cut==5){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    savepath += "dR_nbjets/";
    ischi2 = 0;
  }
  
  Logger logger_class(savepath);
  logger_class.start_log(sample, isPU, ischi2, additional_cut); // write logfile
  
  /// [0] = Signal ///
  chain_mu[0] = new TChain("Muon");
  chain_event[0] = new TChain("Event");
  chain_jet[0] = new TChain("Jet");
  chain_mu[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
  chain_event[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
  chain_jet[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
 
  /// [1] = Data ///
  chain_mu[1] = new TChain("Muon");
  chain_event[1] = new TChain("Event");
  chain_jet[1] = new TChain("Jet");
  chain_mu[1]->Add("./files/data/*.root");
  chain_event[1]->Add("./files/data/*.root");
  chain_jet[1]->Add("./files/data/*.root");
  
  /// [2]~ = Bkg ///
  TString chainlist[N_chain] = {"signal", "data",
    "DY10to50", "DY50plus", "Zbb",  // DY
    "WZtollln_mg", "WZtollqq_mg", "WZtoqqln_mg", "ZZtollll_mg", "ZZtollnn_mg", "ZZtollqq_mg", "WW_mg",  // VV
    "Wbb", "topDIL", "TTG", "TTWW", "WWG", "WWW", "WWZ", "WZZ", "ZZZ", "ttZ", // others (first twos are Wjetes and top)
    "HtoWW", "ggHtoZZ", // Higgs
    "Wtollln_new"}; // Wtollln
  for(int i=2; i<N_chain; i++){
    chain_mu[i] = new TChain("Muon");
    chain_event[i] = new TChain("Event");
    chain_jet[i] = new TChain("Jet");
    chain_mu[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
    chain_event[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
    chain_jet[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
  }
  
  Muon *mu[N_chain];
  Event *evt[N_chain];
  Jet *jet[N_chain];

  double MET, METphi, nu_Pz, W_mass_fit, MET_fit, weight, PUweight_sum;
  int N_entry, solution_selection, smaller;
  struct Fit fitresult[2];
  Color_t fillcolors[N_chain] = {0, kBlack, // signal and data
    kAzure+8, kAzure+8, kAzure+8, // DY
    kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, // VV
    kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, // others
    kViolet, kViolet, // higgs
    kYellow}; // Wtollln
  
  TH1D *hist_met[N_chain];
  TH1D *hist_lead_pt[N_chain];
  TH1D *hist_lead_eta[N_chain];
  TH1D *hist_n_jets[N_chain];
  TH1D *hist_n_bjets[N_chain];
  TH1D *hist_met_err = new TH1D("hist_met_err", "", (met_x_max-met_x_min)/met_dx, met_x_min, met_x_max);
  TH1D *hist_lead_pt_err = new TH1D("hist_lead_pt_err", "", (lead_pt_x_max-lead_pt_x_min)/lead_pt_dx, lead_pt_x_min, lead_pt_x_max);
  TH1D *hist_lead_eta_err = new TH1D("hist_lead_eta_err", "", (lead_eta_x_max-lead_eta_x_min)/lead_eta_dx, lead_eta_x_min, lead_eta_x_max);
  TH1D *hist_n_jets_err = new TH1D("hist_n_jets_err", "", (n_jets_x_max-n_jets_x_min)/n_jets_dx, n_jets_x_min, n_jets_x_max);
  TH1D *hist_n_bjets_err = new TH1D("hist_n_bjets_err", "", (n_bjets_x_max-n_bjets_x_min)/n_bjets_dx, n_bjets_x_min, n_bjets_x_max);
  THStack *bkgstack_met = new THStack("bkgstack_met", "");
  THStack *bkgstack_lead_pt = new THStack("bkgstack_lead_pt", "");
  THStack *bkgstack_lead_eta = new THStack("bkgstack_lead_eta", "");
  THStack *bkgstack_n_jets = new THStack("bkgstack_n_jets", "");
  THStack *bkgstack_n_bjets = new THStack("bkgstack_n_bjets", "");
  
  TH1D *hist_sig = new TH1D("N_sig", "N_sig", 1, 0, 1);
  TH1D *hist_dat = new TH1D("N_dat", "N_dat", 1, 0, 1);
  TH1D *hist_bkg = new TH1D("N_bkg", "N_bkg", 1, 0, 1);
  
  TLorentzVector lep[3], W_real, HN, nu, gen_nu, gamma_star, z_candidate, selection_nu[2];

  int FillOrder[N_chain] = { 0, 1, // signal, data
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, // others
    22, 23, // Higgs
    24, // Wtollln
    5, 6, 7, 8, 9, 10, 11, // VV
    2, 3, 4 }; // DY
  
  for(int l=0; l<N_chain; l++){
    int k = FillOrder[l];
    N_entry = 0;
    PUweight_sum = 0;
 
    mu[k] = new Muon(chain_mu[k], 3);
    evt[k] = new Event(chain_event[k]);
    jet[k] = new Jet(chain_jet[k], 10);
    
    hist_met[k] = new TH1D("hist_met"+TString::Itoa(k,10), "", (met_x_max-met_x_min)/met_dx, met_x_min, met_x_max);
    hist_lead_pt[k] = new TH1D("hist_lead_pt"+TString::Itoa(k,10), "", (lead_pt_x_max-lead_pt_x_min)/lead_pt_dx, lead_pt_x_min, lead_pt_x_max);
    hist_lead_eta[k] = new TH1D("hist_lead_eta"+TString::Itoa(k,10), "", (lead_eta_x_max-lead_eta_x_min)/lead_eta_dx, lead_eta_x_min, lead_eta_x_max);
    hist_n_jets[k] = new TH1D("hist_n_jets"+TString::Itoa(k,10), "", (n_jets_x_max-n_jets_x_min)/n_jets_dx, n_jets_x_min, n_jets_x_max);
    hist_n_bjets[k] = new TH1D("hist_n_bjets"+TString::Itoa(k,10), "", (n_bjets_x_max-n_bjets_x_min)/n_bjets_dx, n_bjets_x_min, n_bjets_x_max);

    cout << endl << "Running Sample : " << chainlist[k] << endl;
    
    for(int i=0; i<chain_mu[k]->GetEntries(); i++){
      
      printeventnumber(i, chain_mu[k]->GetEntries());
      
      chain_mu[k]->GetEntry(i);
      chain_event[k]->GetEntry(i);
      chain_jet[k]->GetEntry(i);
      
      for(int j=0; j<3; j++){ // lepton 4-vec
        lep[j].SetPxPyPzE(mu[k]->Px[j],mu[k]->Py[j],mu[k]->Pz[j],mu[k]->E[j]);
      }
      
      //Cut : Three HNLoose//
      if(muonidsel == "loose")  selection = mu[k]->isHNLooseMuon(0, 15) && mu[k]->isHNLooseMuon(1, 10) && mu[k]->isHNLooseMuon(2, 10);
      if(muonidsel == "tight")  selection = mu[k]->isHNTightMuon(0, 15) && mu[k]->isHNTightMuon(1, 10) && mu[k]->isHNTightMuon(2, 10);
      //Cut : dR > 0.5//
      if(additional_cut == 1 || additional_cut == 2 || additional_cut == 3 || additional_cut == 4 || additional_cut == 5){
        selection *= TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ) > 0.5;
      }
      
      if( selection ){
        MET = evt[k]->MET;
        METphi = evt[k]->METphi;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        
        if(ischi2){ // chi2 fitting
          fit(&fitresult[0], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "m");
          fit(&fitresult[1], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 0 = minus, 1 = plus
          PutNuPz(&selection_nu[0], fitresult[0].Pz_fit);
          PutNuPz(&selection_nu[1], fitresult[1].Pz_fit);
        }
        else{
          PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
          PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus
        }
        
        if( selection_nu[0].Pz() == selection_nu[1].Pz() ){ // solution selection => smaller
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
                
        PutNuPz(&nu, selection_nu[solution_selection].Pz()); // reconstruct HN and W_real 4-vec with selected Pz solution
        if(sample == 60){
          HN = lep[1] + lep[0] + nu;
          W_real = HN + lep[2];
        }
        else{
          HN = lep[1] + lep[2] + nu;
          W_real = HN + lep[0];
        }
        
        //Cut : W(Reco) < 100 GeV//
        if(additional_cut == 2){
          selection *= W_real.M() < 100.0;
        }
        //Cut : chi2 < 0.1//
        if(additional_cut == 3 || additional_cut == 3){
          selection *= fitresult[solution_selection].chi2 < 0.1;
        }
        //Cut : lambda = 0//
        if(additional_cut == 4){
          selection *= fitresult[solution_selection].lambda_fit == 0;
        }
        //Cut : n_bjetes >=1//
        if(additional_cut == 5){
          selection *= evt[k]->n_bjet >= 1;
        }
        
        if(selection){
          N_entry++; // pass-cut calculation
          
          // gamma_star, z_candidate 4-vec reconstruction //
          if( lep[0].DeltaR(lep[1]) < lep[2].DeltaR(lep[1]) ) gamma_star = lep[1] + lep[0];
          else                                                gamma_star = lep[1] + lep[2];
          if( fabs( (lep[1] + lep[0]).M()-91.19 ) >= fabs( (lep[1] + lep[2]).M()-91.19 ) ) z_candidate = lep[1] + lep[2];
          else                                                                             z_candidate = lep[1] + lep[0];
         
          // weight, PU total effect //
          if(!isPU) weight = evt[k]->weight;
          else      weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
          PUweight_sum+=evt[k]->PU_reweight;
          
          // histogram filling //
          int lead_lep_index = findleadlep(lep);
          //if(z_candidate.M() > 50 && z_candidate.M() < 70){
          hist_met[k]->Fill(MET, weight);
          //}
          hist_lead_pt[k]->Fill(lep[lead_lep_index].Pt(), weight);
          hist_lead_eta[k]->Fill(lep[lead_lep_index].Eta(), weight);
          hist_n_jets[k]->Fill(jet[k]->n_Jet, weight);
          hist_n_bjets[k]->Fill(evt[k]->n_bjet, weight);
          
          // one bin filling //
          if(k==0){
            hist_sig->Fill(0., weight);
          }
          else if(k==1){
            hist_dat->Fill(0., weight);
          }
          else{
            hist_bkg->Fill(0., weight);
          }
          
        } // additional cut loop
      } // HNLoose + deltaR cut loop
    } // event loop

    cout << "weight = " << evt[k]->weight << endl
         << "entry = " << N_entry << endl
    << "PUweight = " << PUweight_sum/(double)N_entry << endl;
    
    logger_class.write_n_events(chainlist[k], evt[k]->weight, N_entry, PUweight_sum/(double)N_entry);
    
    hist_met[k]->SetLineColor(fillcolors[k]);
    hist_lead_pt[k]->SetLineColor(fillcolors[k]);
    hist_lead_eta[k]->SetLineColor(fillcolors[k]);
    hist_n_jets[k]->SetLineColor(fillcolors[k]);
    hist_n_bjets[k]->SetLineColor(fillcolors[k]);
    hist_met[k]->SetFillColor(fillcolors[k]);
    hist_lead_pt[k]->SetFillColor(fillcolors[k]);
    hist_lead_eta[k]->SetFillColor(fillcolors[k]);
    hist_n_jets[k]->SetFillColor(fillcolors[k]);
    hist_n_bjets[k]->SetFillColor(fillcolors[k]);
    
    /// make TH1 for error & make THSatck
    if(k>=2){
      bkgstack_met->Add(hist_met[k]);
      bkgstack_lead_pt->Add(hist_lead_pt[k]);
      bkgstack_lead_eta->Add(hist_lead_eta[k]);
      bkgstack_n_jets->Add(hist_n_jets[k]);
      bkgstack_n_bjets->Add(hist_n_bjets[k]);
      hist_met_err->Add(hist_met[k]);
      hist_lead_pt_err->Add(hist_lead_pt[k]);
      hist_lead_eta_err->Add(hist_lead_eta[k]);
      hist_n_jets_err->Add(hist_n_jets[k]);
      hist_n_bjets_err->Add(hist_n_bjets[k]);
    }

  
  } // chain loop
  
  logger_class.divide_log();
  
  /// signal scale factor ///
  double k_factor = 1.34, coupling_const;
  if(sample==40){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==50){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==60){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==100){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=1;
    }
    else{
      coupling_const=0.01;
    }
  }
  hist_met[0]->Scale(k_factor*coupling_const);
  hist_lead_pt[0]->Scale(k_factor*coupling_const);
  hist_lead_eta[0]->Scale(k_factor*coupling_const);
  hist_n_jets[0]->Scale(k_factor*coupling_const);
  hist_n_bjets[0]->Scale(k_factor*coupling_const);
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  //legend//
  TLegend *lg1 = new TLegend(0.65, 0.55, 0.97, 0.90);
  lg1->SetFillStyle(0);
  lg1->SetBorderSize(0);
  lg1->AddEntry(hist_met[1], "data", "p");
  lg1->AddEntry(hist_met[2], "DY", "f");
  lg1->AddEntry(hist_met[5], "VV", "f");
  lg1->AddEntry(hist_met[24], "W#rightarrowlll#nu", "f");
  lg1->AddEntry(hist_met[22], "Higgs", "f");
  lg1->AddEntry(hist_met[14], "others", "f");
  lg1->AddEntry(hist_met[0], "HN"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  //y=0//
  double x_0[2], y_0[2];
  x_0[0] = -1000;  y_0[0] = 0;
  x_0[1] = 1000;  y_0[1] = 0;
  TGraph *g0 = new TGraph(2, x_0, y_0);
  //y=1//
  double x_1[2], y_1[2];
  x_1[0] = -1000;  y_1[0] = 1;
  x_1[1] = 1000;  y_1[1] = 1;
  TGraph *g1 = new TGraph(2, x_1, y_1);
  
  /// met ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  TPad *c1_up = new TPad("", "", 0, 0.25, 1, 1);
  c1_up->SetTopMargin( 0.05 ); c1_up->SetBottomMargin( 0.02 ); c1_up->SetRightMargin( 0.02 ); c1_up->SetLeftMargin( 0.1 );
  TPad *c1_down = new TPad("", "", 0, 0, 1, 0.25);
  c1_down->SetTopMargin( 0.03 ); c1_down->SetBottomMargin( 0.25 ); c1_down->SetRightMargin( 0.02 ); c1_down->SetLeftMargin( 0.1 ); c1_down->SetGridx(); c1_down->SetGridy();
  c1_up->Draw();
  c1_down->Draw();
  c1_up->cd();
  ///bkg///
  bkgstack_met->Draw("hist");
  bkgstack_met->SetMaximum(y_max[0]); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_met->GetXaxis()->SetLabelSize(0);
  bkgstack_met->GetYaxis()->SetTitle("Events / "+TString::Itoa(met_dx,10)+" GeV");
  bkgstack_met->GetYaxis()->SetLabelSize(0.04);
  bkgstack_met->GetYaxis()->SetTitleSize(0.04);
  bkgstack_met->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_met[1]->SetMarkerStyle(3);
  hist_met[1]->SetMarkerSize(1);
  hist_met[1]->Draw("APsameE1");
  ///signal///
  hist_met[0]->SetLineColor(kRed);
  hist_met[0]->Draw("histsame");
  ///err///
  hist_met_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_met_err->SetFillStyle(3004);
  hist_met_err->SetFillColor(kBlue);
  hist_met_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c1_down->cd();
  TH1D *hist_met_diff = (TH1D*)hist_met[1]->Clone();
  hist_met_diff->Divide(hist_met_err);
  //hist_met_diff->Add(hist_met_err,-1);
  hist_met_diff->SetMaximum(2);
  hist_met_diff->SetMinimum(0);
  hist_met_diff->GetXaxis()->SetTitle("MET [GeV]");
  hist_met_diff->GetXaxis()->SetLabelSize(0.10);
  hist_met_diff->GetXaxis()->SetTitleSize(0.10);
  hist_met_diff->GetYaxis()->SetLabelSize(0.08);
  hist_met_diff->SetYTitle("#frac{DATA}{MC}");
  hist_met_diff->GetYaxis()->SetTitleSize(0.08);
  hist_met_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_met_diff->SetFillColorAlpha(45,0.35);
  hist_met_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");

  /// lead_pt ///
  TCanvas *c2 = new TCanvas("c2", "", 800, 800);
  TPad *c2_up = new TPad("", "", 0, 0.25, 1, 1);
  c2_up->SetTopMargin( 0.05 ); c2_up->SetBottomMargin( 0.02 ); c2_up->SetRightMargin( 0.02 ); c2_up->SetLeftMargin( 0.1 );
  TPad *c2_down = new TPad("", "", 0, 0, 1, 0.25);
  c2_down->SetTopMargin( 0.03 ); c2_down->SetBottomMargin( 0.25 ); c2_down->SetRightMargin( 0.02 ); c2_down->SetLeftMargin( 0.1 ); c2_down->SetGridx(); c2_down->SetGridy();
  c2_up->Draw();
  c2_down->Draw();
  c2_up->cd();
  ///bkg///
  bkgstack_lead_pt->Draw("hist");
  bkgstack_lead_pt->SetMaximum(y_max[1]);
  bkgstack_lead_pt->GetXaxis()->SetLabelSize(0);
  bkgstack_lead_pt->GetYaxis()->SetTitle("Events / "+TString::Itoa(lead_pt_dx,10)+" GeV");
  bkgstack_lead_pt->GetYaxis()->SetLabelSize(0.04);
  bkgstack_lead_pt->GetYaxis()->SetTitleSize(0.04);
  bkgstack_lead_pt->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_lead_pt[1]->SetMarkerStyle(3);
  hist_lead_pt[1]->SetMarkerSize(1);
  hist_lead_pt[1]->Draw("APsameE1");
  ///signal///
  hist_lead_pt[0]->SetLineColor(kRed);
  hist_lead_pt[0]->Draw("histsame");
  ///err///
  hist_lead_pt_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_lead_pt_err->SetFillStyle(3004);
  hist_lead_pt_err->SetFillColor(kBlue);
  hist_lead_pt_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c2_down->cd();
  TH1D *hist_lead_pt_diff = (TH1D*)hist_lead_pt[1]->Clone();
  hist_lead_pt_diff->Divide(hist_lead_pt_err);
  //hist_lead_pt_diff->Add(hist_lead_pt_err,-1);
  hist_lead_pt_diff->SetMaximum(2);
  hist_lead_pt_diff->SetMinimum(0);
  hist_lead_pt_diff->GetXaxis()->SetTitle("P_{T} of the leading lepton [GeV]");
  hist_lead_pt_diff->GetXaxis()->SetLabelSize(0.10);
  hist_lead_pt_diff->GetXaxis()->SetTitleSize(0.10);
  hist_lead_pt_diff->GetYaxis()->SetLabelSize(0.08);
  hist_lead_pt_diff->SetYTitle("#frac{DATA}{MC}");
  hist_lead_pt_diff->GetYaxis()->SetTitleSize(0.08);
  hist_lead_pt_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_lead_pt_diff->SetFillColorAlpha(45,0.35);
  hist_lead_pt_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// lead_eta ///
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  TPad *c3_up = new TPad("", "", 0, 0.25, 1, 1);
  c3_up->SetTopMargin( 0.05 ); c3_up->SetBottomMargin( 0.02 ); c3_up->SetRightMargin( 0.02 ); c3_up->SetLeftMargin( 0.1 );
  TPad *c3_down = new TPad("", "", 0, 0, 1, 0.25);
  c3_down->SetTopMargin( 0.03 ); c3_down->SetBottomMargin( 0.25 ); c3_down->SetRightMargin( 0.02 ); c3_down->SetLeftMargin( 0.1 ); c3_down->SetGridx(); c3_down->SetGridy();
  c3_up->Draw();
  c3_down->Draw();
  c3_up->cd();
  ///bkg///
  bkgstack_lead_eta->Draw("hist");
  bkgstack_lead_eta->SetMaximum(y_max[2]);
  bkgstack_lead_eta->GetXaxis()->SetLabelSize(0);
  bkgstack_lead_eta->GetYaxis()->SetTitle("Events / 0."+TString::Itoa(10*lead_eta_dx,10));
  bkgstack_lead_eta->GetYaxis()->SetLabelSize(0.04);
  bkgstack_lead_eta->GetYaxis()->SetTitleSize(0.04);
  bkgstack_lead_eta->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_lead_eta[1]->SetMarkerStyle(3);
  hist_lead_eta[1]->SetMarkerSize(1);
  hist_lead_eta[1]->Draw("APsameE1");
  ///signal///
  hist_lead_eta[0]->SetLineColor(kRed);
  hist_lead_eta[0]->Draw("histsame");
  ///err///
  hist_lead_eta_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_lead_eta_err->SetFillStyle(3004);
  hist_lead_eta_err->SetFillColor(kBlue);
  hist_lead_eta_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c3_down->cd();
  TH1D *hist_lead_eta_diff = (TH1D*)hist_lead_eta[1]->Clone();
  hist_lead_eta_diff->Divide(hist_lead_eta_err);
  //hist_lead_eta_diff->Add(hist_lead_eta_err,-1);
  hist_lead_eta_diff->SetMaximum(2);
  hist_lead_eta_diff->SetMinimum(0);
  hist_lead_eta_diff->GetXaxis()->SetTitle("#eta of the leading lepton");
  hist_lead_eta_diff->GetXaxis()->SetLabelSize(0.10);
  hist_lead_eta_diff->GetXaxis()->SetTitleSize(0.10);
  hist_lead_eta_diff->GetYaxis()->SetLabelSize(0.08);
  hist_lead_eta_diff->SetYTitle("#frac{DATA}{MC}");
  hist_lead_eta_diff->GetYaxis()->SetTitleSize(0.08);
  hist_lead_eta_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_lead_eta_diff->SetFillColorAlpha(45,0.35);
  hist_lead_eta_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// n_jets ///
  TCanvas *c4 = new TCanvas("c4", "", 800, 800);
  TPad *c4_up = new TPad("", "", 0, 0.25, 1, 1);
  c4_up->SetTopMargin( 0.05 ); c4_up->SetBottomMargin( 0.02 ); c4_up->SetRightMargin( 0.02 ); c4_up->SetLeftMargin( 0.1 );
  TPad *c4_down = new TPad("", "", 0, 0, 1, 0.25);
  c4_down->SetTopMargin( 0.03 ); c4_down->SetBottomMargin( 0.25 ); c4_down->SetRightMargin( 0.02 ); c4_down->SetLeftMargin( 0.1 ); c4_down->SetGridx(); c4_down->SetGridy();
  c4_up->Draw();
  c4_down->Draw();
  c4_up->cd();
  ///bkg///
  bkgstack_n_jets->Draw("hist");
  bkgstack_n_jets->SetMaximum(y_max[3]);
  bkgstack_n_jets->GetXaxis()->SetLabelSize(0);
  bkgstack_n_jets->GetYaxis()->SetTitle("Events");
  bkgstack_n_jets->GetYaxis()->SetLabelSize(0.04);
  bkgstack_n_jets->GetYaxis()->SetTitleSize(0.04);
  bkgstack_n_jets->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_n_jets[1]->SetMarkerStyle(3);
  hist_n_jets[1]->SetMarkerSize(1);
  hist_n_jets[1]->Draw("APsameE1");
  ///signal///
  hist_n_jets[0]->SetLineColor(kRed);
  hist_n_jets[0]->Draw("histsame");
  ///err///
  hist_n_jets_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_n_jets_err->SetFillStyle(3004);
  hist_n_jets_err->SetFillColor(kBlue);
  hist_n_jets_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c4_down->cd();
  TH1D *hist_n_jets_diff = (TH1D*)hist_n_jets[1]->Clone();
  hist_n_jets_diff->Divide(hist_n_jets_err);
  //hist_n_jets_diff->Add(hist_n_jets_err,-1);
  hist_n_jets_diff->SetMaximum(2);
  hist_n_jets_diff->SetMinimum(0);
  hist_n_jets_diff->GetXaxis()->SetTitle("# of jets");
  hist_n_jets_diff->GetXaxis()->SetLabelSize(0.10);
  hist_n_jets_diff->GetXaxis()->SetTitleSize(0.10);
  hist_n_jets_diff->GetYaxis()->SetLabelSize(0.08);
  hist_n_jets_diff->SetYTitle("#frac{DATA}{MC}");
  hist_n_jets_diff->GetYaxis()->SetTitleSize(0.08);
  hist_n_jets_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_n_jets_diff->SetFillColorAlpha(45,0.35);
  hist_n_jets_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// n_bjets ///
  TCanvas *c5 = new TCanvas("c5", "", 800, 800);
  TPad *c5_up = new TPad("", "", 0, 0.25, 1, 1);
  c5_up->SetTopMargin( 0.05 ); c5_up->SetBottomMargin( 0.02 ); c5_up->SetRightMargin( 0.02 ); c5_up->SetLeftMargin( 0.1 );
  TPad *c5_down = new TPad("", "", 0, 0, 1, 0.25);
  c5_down->SetTopMargin( 0.03 ); c5_down->SetBottomMargin( 0.25 ); c5_down->SetRightMargin( 0.02 ); c5_down->SetLeftMargin( 0.1 ); c5_down->SetGridx(); c5_down->SetGridy();
  c5_up->Draw();
  c5_down->Draw();
  c5_up->cd();
  c5_up->SetLogy();
  ///bkg///
  bkgstack_n_bjets->Draw("hist");
  bkgstack_n_bjets->SetMaximum(y_max[4]);
  bkgstack_n_bjets->GetXaxis()->SetLabelSize(0);
  bkgstack_n_bjets->GetYaxis()->SetTitle("Events");
  bkgstack_n_bjets->GetYaxis()->SetLabelSize(0.04);
  bkgstack_n_bjets->GetYaxis()->SetTitleSize(0.04);
  bkgstack_n_bjets->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_n_bjets[1]->SetMarkerStyle(3);
  hist_n_bjets[1]->SetMarkerSize(1);
  hist_n_bjets[1]->Draw("APsameE1");
  ///signal///
  hist_n_bjets[0]->SetLineColor(kRed);
  hist_n_bjets[0]->Draw("histsame");
  ///err///
  hist_n_bjets_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_n_bjets_err->SetFillStyle(3004);
  hist_n_bjets_err->SetFillColor(kBlue);
  hist_n_bjets_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c5_down->cd();
  TH1D *hist_n_bjets_diff = (TH1D*)hist_n_bjets[1]->Clone();
  hist_n_bjets_diff->Divide(hist_n_bjets_err);
  //hist_n_bjets_diff->Add(hist_n_bjets_err,-1);
  hist_n_bjets_diff->SetMaximum(2);
  hist_n_bjets_diff->SetMinimum(0);
  hist_n_bjets_diff->GetXaxis()->SetTitle("# of b-tagged jets");
  hist_n_bjets_diff->GetXaxis()->SetLabelSize(0.10);
  hist_n_bjets_diff->GetXaxis()->SetTitleSize(0.10);
  hist_n_bjets_diff->GetYaxis()->SetLabelSize(0.08);
  hist_n_bjets_diff->SetYTitle("#frac{DATA}{MC}");
  hist_n_bjets_diff->GetYaxis()->SetTitleSize(0.08);
  hist_n_bjets_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_n_bjets_diff->SetFillColorAlpha(45,0.35);
  hist_n_bjets_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  if(isSave){
    c1->SaveAs(savepath+TString::Itoa(sample, 10)+"_met.root");
    c2->SaveAs(savepath+TString::Itoa(sample, 10)+"_lead_pt.root");
    c3->SaveAs(savepath+TString::Itoa(sample, 10)+"_lead_eta.root");
    c4->SaveAs(savepath+TString::Itoa(sample, 10)+"_n_jets.root");
    c5->SaveAs(savepath+TString::Itoa(sample, 10)+"_n_bjets.root");
  }
  if(isClose){
    c1->Close();
    c2->Close();
    c3->Close();
    c4->Close();
    c5->Close();
  }
  cout
  << endl
  << "Logfile ------> " << logger_class.GetLogfilePath() << endl;
  logger_class.write_onebin(hist_sig, hist_dat, hist_bkg);
}

//////////////////////////////////
/// main function (for script) ///
//////////////////////////////////

void control_check_mumumu(int sample, int isPU, int additional_cut, int muonidsel_int){ // muonidsel = 0 = loose, muonidsel = 1 = tight
  
  TH1::SetDefaultSumw2(true);
  
  TStopwatch totaltime;
  totaltime.Start();
  //int sample, isPU, ischi2, additional_cut
  int isSave=1, isClose=1, y_max[N_plots], ischi2;
  double met_x_min, met_x_max, met_dx, lead_pt_x_min, lead_pt_x_max, lead_pt_dx, lead_eta_x_min, lead_eta_x_max, lead_eta_dx, n_jets_x_min, n_jets_x_max, n_jets_dx, n_bjets_x_min, n_bjets_x_max, n_bjets_dx;
  TString savepath = "/home/jskim/Documents/CMS/trilepton/mumumu/plots/control/", muonidsel;
  bool selection;
  TChain *chain_mu[N_chain], *chain_event[N_chain], *chain_jet[N_chain];
  /*
   cout
   << "===============================================" << endl
   << "40, 50, 60, 100(no gen lvl info) : signal" << endl
   << "===============================================" << endl
   << "Run : ";
   cin >> sample;
   cout << endl;
   
   cout
   << "===============================================" << endl
   << "PileUP? (0 = no, 1 = yes) : ";
   cin >> isPU;
   cout
   << "chi2 fit? (0 = no, 1 = yes) :";
   cin >> ischi2;
   cout
   << "===============================================" << endl
   << "Additional Cut?" << endl
   << "  0 : No Additional Cut" << endl
   << "  1 : delta R > 0.5" << endl
   << "  2 : delta R > 0.5 & W_reco < 100 GeV" << endl
   << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
   << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda = 0" << endl
   << "  5 : delta R > 0.5 & # of b-tagged jet >= 1" << endl
   << "===============================================" << endl
   << "Cut : ";
   cin >> additional_cut;
   cout << endl;
   */
  
  if(muonidsel_int == 0)  muonidsel = "loose";
  if(muonidsel_int == 1)  muonidsel = "tight";
  
  if(additional_cut==0){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=1000;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=1000;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=1000;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "before_dR/loose_";
    if(muonidsel == "tight")  savepath += "before_dR/tight_";
    ischi2 = 0;
  }
  if(additional_cut==1){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=800;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=800;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=1000;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "dR/loose_";
    if(muonidsel == "tight")  savepath += "dR/tight_";
    ischi2 = 0;
  }
  if(additional_cut==2){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=100;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "dR_W/loose_";
    if(muonidsel == "tight")  savepath += "dR_W/tight_";
    ischi2 = 0;
  }
  if(additional_cut==3){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "dR_chi2/loose_";
    if(muonidsel == "tight")  savepath += "dR_chi2/tight_";
    ischi2 = 1;
  }
  if(additional_cut==4){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "dR_chi2_lambda/loose_";
    if(muonidsel == "tight")  savepath += "dR_chi2_lambda/tight_";
    ischi2 = 1;
  }
  if(additional_cut==5){
    met_x_min=0; met_x_max=150; met_dx=10; y_max[0]=60;
    lead_pt_x_min=0; lead_pt_x_max=200; lead_pt_dx=10; y_max[1]=100;
    lead_eta_x_min=-2.5; lead_eta_x_max=2.5; lead_eta_dx=0.5; y_max[2]=100;
    n_jets_x_min=0; n_jets_x_max=10; n_jets_dx=1; y_max[3]=1500;
    n_bjets_x_min=0; n_bjets_x_max=5; n_bjets_dx=1; y_max[4]=3000;
    if(muonidsel == "loose")  savepath += "dR_nbjets/loose_";
    if(muonidsel == "tight")  savepath += "dR_nbjets/tight_";
    ischi2 = 0;
  }
  
  Logger logger_class(savepath);
  logger_class.start_log(sample, isPU, ischi2, additional_cut); // write logfile
  
  /// [0] = Signal ///
  chain_mu[0] = new TChain("Muon");
  chain_event[0] = new TChain("Event");
  chain_jet[0] = new TChain("Jet");
  chain_mu[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
  chain_event[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
  chain_jet[0]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample, 10)+"_new_5_3_14.root");
  
  /// [1] = Data ///
  chain_mu[1] = new TChain("Muon");
  chain_event[1] = new TChain("Event");
  chain_jet[1] = new TChain("Jet");
  chain_mu[1]->Add("./files/data/*.root");
  chain_event[1]->Add("./files/data/*.root");
  chain_jet[1]->Add("./files/data/*.root");
  
  /// [2]~ = Bkg ///
  TString chainlist[N_chain] = {"signal", "data",
    "DY10to50", "DY50plus", "Zbb",  // DY
    "WZtollln_mg", "WZtollqq_mg", "WZtoqqln_mg", "ZZtollll_mg", "ZZtollnn_mg", "ZZtollqq_mg", "WW_mg",  // VV
    "Wbb", "topDIL", "TTG", "TTWW", "WWG", "WWW", "WWZ", "WZZ", "ZZZ", "ttZ", // others (first twos are Wjetes and top)
    "HtoWW", "ggHtoZZ", // Higgs
    "Wtollln_new"}; // Wtollln
  for(int i=2; i<N_chain; i++){
    chain_mu[i] = new TChain("Muon");
    chain_event[i] = new TChain("Event");
    chain_jet[i] = new TChain("Jet");
    chain_mu[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
    chain_event[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
    chain_jet[i]->Add("./files/bkg/trilepton_mumumu_SK"+chainlist[i]+"_5_3_14.root");
  }
  
  Muon *mu[N_chain];
  Event *evt[N_chain];
  Jet *jet[N_chain];
  
  double MET, METphi, nu_Pz, W_mass_fit, MET_fit, weight, PUweight_sum;
  int N_entry, solution_selection, smaller;
  struct Fit fitresult[2];
  Color_t fillcolors[N_chain] = {0, kBlack, // signal and data
    kAzure+8, kAzure+8, kAzure+8, // DY
    kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, // VV
    kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, // others
    kViolet, kViolet, // higgs
    kYellow}; // Wtollln
  
  TH1D *hist_met[N_chain];
  TH1D *hist_lead_pt[N_chain];
  TH1D *hist_lead_eta[N_chain];
  TH1D *hist_n_jets[N_chain];
  TH1D *hist_n_bjets[N_chain];
  TH1D *hist_met_err = new TH1D("hist_met_err", "", (met_x_max-met_x_min)/met_dx, met_x_min, met_x_max);
  TH1D *hist_lead_pt_err = new TH1D("hist_lead_pt_err", "", (lead_pt_x_max-lead_pt_x_min)/lead_pt_dx, lead_pt_x_min, lead_pt_x_max);
  TH1D *hist_lead_eta_err = new TH1D("hist_lead_eta_err", "", (lead_eta_x_max-lead_eta_x_min)/lead_eta_dx, lead_eta_x_min, lead_eta_x_max);
  TH1D *hist_n_jets_err = new TH1D("hist_n_jets_err", "", (n_jets_x_max-n_jets_x_min)/n_jets_dx, n_jets_x_min, n_jets_x_max);
  TH1D *hist_n_bjets_err = new TH1D("hist_n_bjets_err", "", (n_bjets_x_max-n_bjets_x_min)/n_bjets_dx, n_bjets_x_min, n_bjets_x_max);
  THStack *bkgstack_met = new THStack("bkgstack_met", "");
  THStack *bkgstack_lead_pt = new THStack("bkgstack_lead_pt", "");
  THStack *bkgstack_lead_eta = new THStack("bkgstack_lead_eta", "");
  THStack *bkgstack_n_jets = new THStack("bkgstack_n_jets", "");
  THStack *bkgstack_n_bjets = new THStack("bkgstack_n_bjets", "");
  
  TH1D *hist_sig = new TH1D("N_sig", "N_sig", 1, 0, 1);
  TH1D *hist_dat = new TH1D("N_dat", "N_dat", 1, 0, 1);
  TH1D *hist_bkg = new TH1D("N_bkg", "N_bkg", 1, 0, 1);
  
  TLorentzVector lep[3], W_real, HN, nu, gen_nu, gamma_star, z_candidate, selection_nu[2];
  
  int FillOrder[N_chain] = { 0, 1, // signal, data
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, // others
    22, 23, // Higgs
    24, // Wtollln
    5, 6, 7, 8, 9, 10, 11, // VV
    2, 3, 4 }; // DY
  
  for(int l=0; l<N_chain; l++){
    int k = FillOrder[l];
    N_entry = 0;
    PUweight_sum = 0;
    
    mu[k] = new Muon(chain_mu[k], 3);
    evt[k] = new Event(chain_event[k]);
    jet[k] = new Jet(chain_jet[k], 10);
    
    hist_met[k] = new TH1D("hist_met"+TString::Itoa(k,10), "", (met_x_max-met_x_min)/met_dx, met_x_min, met_x_max);
    hist_lead_pt[k] = new TH1D("hist_lead_pt"+TString::Itoa(k,10), "", (lead_pt_x_max-lead_pt_x_min)/lead_pt_dx, lead_pt_x_min, lead_pt_x_max);
    hist_lead_eta[k] = new TH1D("hist_lead_eta"+TString::Itoa(k,10), "", (lead_eta_x_max-lead_eta_x_min)/lead_eta_dx, lead_eta_x_min, lead_eta_x_max);
    hist_n_jets[k] = new TH1D("hist_n_jets"+TString::Itoa(k,10), "", (n_jets_x_max-n_jets_x_min)/n_jets_dx, n_jets_x_min, n_jets_x_max);
    hist_n_bjets[k] = new TH1D("hist_n_bjets"+TString::Itoa(k,10), "", (n_bjets_x_max-n_bjets_x_min)/n_bjets_dx, n_bjets_x_min, n_bjets_x_max);
    
    cout << endl << "Running Sample : " << chainlist[k] << endl;
    
    for(int i=0; i<chain_mu[k]->GetEntries(); i++){
      
      printeventnumber(i, chain_mu[k]->GetEntries());
      
      chain_mu[k]->GetEntry(i);
      chain_event[k]->GetEntry(i);
      chain_jet[k]->GetEntry(i);
      
      for(int j=0; j<3; j++){ // lepton 4-vec
        lep[j].SetPxPyPzE(mu[k]->Px[j],mu[k]->Py[j],mu[k]->Pz[j],mu[k]->E[j]);
      }
      
      //Cut : Three HNLoose//
      if(muonidsel == "loose")  selection = mu[k]->isHNLooseMuon(0, 15) && mu[k]->isHNLooseMuon(1, 10) && mu[k]->isHNLooseMuon(2, 10);
      if(muonidsel == "tight")  selection = mu[k]->isHNTightMuon(0, 15) && mu[k]->isHNTightMuon(1, 10) && mu[k]->isHNTightMuon(2, 10);
      //Cut : dR > 0.5//
      if(additional_cut == 1 || additional_cut == 2 || additional_cut == 3 || additional_cut == 4 || additional_cut == 5){
        selection *= TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ) > 0.5;
      }
      
      if( selection ){
        MET = evt[k]->MET;
        METphi = evt[k]->METphi;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        
        if(ischi2){ // chi2 fitting
          fit(&fitresult[0], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "m");
          fit(&fitresult[1], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 0 = minus, 1 = plus
          PutNuPz(&selection_nu[0], fitresult[0].Pz_fit);
          PutNuPz(&selection_nu[1], fitresult[1].Pz_fit);
        }
        else{
          PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
          PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus
        }
        
        if( selection_nu[0].Pz() == selection_nu[1].Pz() ){ // solution selection => smaller
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
        
        PutNuPz(&nu, selection_nu[solution_selection].Pz()); // reconstruct HN and W_real 4-vec with selected Pz solution
        if(sample == 60){
          HN = lep[1] + lep[0] + nu;
          W_real = HN + lep[2];
        }
        else{
          HN = lep[1] + lep[2] + nu;
          W_real = HN + lep[0];
        }
        
        //Cut : W(Reco) < 100 GeV//
        if(additional_cut == 2){
          selection *= W_real.M() < 100.0;
        }
        //Cut : chi2 < 0.1//
        if(additional_cut == 3 || additional_cut == 3){
          selection *= fitresult[solution_selection].chi2 < 0.1;
        }
        //Cut : lambda = 0//
        if(additional_cut == 4){
          selection *= fitresult[solution_selection].lambda_fit == 0;
        }
        //Cut : n_bjetes >=1//
        if(additional_cut == 5){
          selection *= evt[k]->n_bjet >= 1;
        }
        
        if(selection){
          N_entry++; // pass-cut calculation
          
          // gamma_star, z_candidate 4-vec reconstruction //
          if( lep[0].DeltaR(lep[1]) < lep[2].DeltaR(lep[1]) ) gamma_star = lep[1] + lep[0];
          else                                                gamma_star = lep[1] + lep[2];
          if( fabs( (lep[1] + lep[0]).M()-91.19 ) >= fabs( (lep[1] + lep[2]).M()-91.19 ) ) z_candidate = lep[1] + lep[2];
          else                                                                             z_candidate = lep[1] + lep[0];
          
          // weight, PU total effect //
          if(!isPU) weight = evt[k]->weight;
          else      weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
          PUweight_sum+=evt[k]->PU_reweight;
          
          // histogram filling //
          int lead_lep_index = findleadlep(lep);
          hist_met[k]->Fill(MET, weight);
          hist_lead_pt[k]->Fill(lep[lead_lep_index].Pt(), weight);
          hist_lead_eta[k]->Fill(lep[lead_lep_index].Eta(), weight);
          hist_n_jets[k]->Fill(jet[k]->n_Jet, weight);
          hist_n_bjets[k]->Fill(evt[k]->n_bjet, weight);
          
          // one bin filling //
          if(k==0){
            hist_sig->Fill(0., weight);
          }
          else if(k==1){
            hist_dat->Fill(0., weight);
          }
          else{
            hist_bkg->Fill(0., weight);
          }
          
        } // additional cut loop
      } // HNLoose + deltaR cut loop
    } // event loop
    
    cout << "weight = " << evt[k]->weight << endl
    << "entry = " << N_entry << endl
    << "PUweight = " << PUweight_sum/(double)N_entry << endl;
    
    logger_class.write_n_events(chainlist[k], evt[k]->weight, N_entry, PUweight_sum/(double)N_entry);
    
    hist_met[k]->SetLineColor(fillcolors[k]);
    hist_lead_pt[k]->SetLineColor(fillcolors[k]);
    hist_lead_eta[k]->SetLineColor(fillcolors[k]);
    hist_n_jets[k]->SetLineColor(fillcolors[k]);
    hist_n_bjets[k]->SetLineColor(fillcolors[k]);
    hist_met[k]->SetFillColor(fillcolors[k]);
    hist_lead_pt[k]->SetFillColor(fillcolors[k]);
    hist_lead_eta[k]->SetFillColor(fillcolors[k]);
    hist_n_jets[k]->SetFillColor(fillcolors[k]);
    hist_n_bjets[k]->SetFillColor(fillcolors[k]);
    
    /// make TH1 for error & make THSatck
    if(k>=2){
      bkgstack_met->Add(hist_met[k]);
      bkgstack_lead_pt->Add(hist_lead_pt[k]);
      bkgstack_lead_eta->Add(hist_lead_eta[k]);
      bkgstack_n_jets->Add(hist_n_jets[k]);
      bkgstack_n_bjets->Add(hist_n_bjets[k]);
      hist_met_err->Add(hist_met[k]);
      hist_lead_pt_err->Add(hist_lead_pt[k]);
      hist_lead_eta_err->Add(hist_lead_eta[k]);
      hist_n_jets_err->Add(hist_n_jets[k]);
      hist_n_bjets_err->Add(hist_n_bjets[k]);
    }
    
    
  } // chain loop
  
  logger_class.divide_log();
  
  /// signal scale factor ///
  double k_factor = 1.34, coupling_const;
  if(sample==40){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==50){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==60){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
  }
  if(sample==100){
    if(additional_cut==0 || additional_cut==1){
      coupling_const=1;
    }
    else{
      coupling_const=0.01;
    }
  }
  hist_met[0]->Scale(k_factor*coupling_const);
  hist_lead_pt[0]->Scale(k_factor*coupling_const);
  hist_lead_eta[0]->Scale(k_factor*coupling_const);
  hist_n_jets[0]->Scale(k_factor*coupling_const);
  hist_n_bjets[0]->Scale(k_factor*coupling_const);
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  //legend//
  TLegend *lg1 = new TLegend(0.65, 0.55, 0.97, 0.90);
  lg1->SetFillStyle(0);
  lg1->SetBorderSize(0);
  lg1->AddEntry(hist_met[1], "data", "p");
  lg1->AddEntry(hist_met[2], "DY", "f");
  lg1->AddEntry(hist_met[5], "VV", "f");
  lg1->AddEntry(hist_met[24], "W#rightarrowlll#nu", "f");
  lg1->AddEntry(hist_met[22], "Higgs", "f");
  lg1->AddEntry(hist_met[14], "others", "f");
  lg1->AddEntry(hist_met[0], "HN"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  //y=0//
  double x_0[2], y_0[2];
  x_0[0] = -1000;  y_0[0] = 0;
  x_0[1] = 1000;  y_0[1] = 0;
  TGraph *g0 = new TGraph(2, x_0, y_0);
  //y=1//
  double x_1[2], y_1[2];
  x_1[0] = -1000;  y_1[0] = 1;
  x_1[1] = 1000;  y_1[1] = 1;
  TGraph *g1 = new TGraph(2, x_1, y_1);
  
  /// met ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  TPad *c1_up = new TPad("c1_up", "", 0, 0.25, 1, 1);
  c1_up->SetTopMargin( 0.05 ); c1_up->SetBottomMargin( 0.02 ); c1_up->SetRightMargin( 0.02 ); c1_up->SetLeftMargin( 0.1 );
  TPad *c1_down = new TPad("c1_down", "", 0, 0, 1, 0.25);
  c1_down->SetTopMargin( 0.03 ); c1_down->SetBottomMargin( 0.25 ); c1_down->SetRightMargin( 0.02 ); c1_down->SetLeftMargin( 0.1 ); c1_down->SetGridx(); c1_down->SetGridy();
  c1_up->Draw();
  c1_down->Draw();
  c1_up->cd();
  ///bkg///
  bkgstack_met->Draw("hist");
  bkgstack_met->SetMaximum(y_max[0]); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_met->GetXaxis()->SetLabelSize(0);
  bkgstack_met->GetYaxis()->SetTitle("Events / "+TString::Itoa(met_dx,10)+" GeV");
  bkgstack_met->GetYaxis()->SetLabelSize(0.04);
  bkgstack_met->GetYaxis()->SetTitleSize(0.04);
  bkgstack_met->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_met[1]->SetMarkerStyle(3);
  hist_met[1]->SetMarkerSize(1);
  hist_met[1]->Draw("APsameE1");
  ///signal///
  hist_met[0]->SetLineColor(kRed);
  hist_met[0]->Draw("histsame");
  ///err///
  hist_met_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_met_err->SetFillStyle(3004);
  hist_met_err->SetFillColor(kBlue);
  hist_met_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c1_down->cd();
  TH1D *hist_met_diff = (TH1D*)hist_met[1]->Clone();
  hist_met_diff->Divide(hist_met_err);
  //hist_met_diff->Add(hist_met_err,-1);
  hist_met_diff->SetMaximum(2);
  hist_met_diff->SetMinimum(0);
  hist_met_diff->GetXaxis()->SetTitle("MET [GeV]");
  hist_met_diff->GetXaxis()->SetLabelSize(0.10);
  hist_met_diff->GetXaxis()->SetTitleSize(0.10);
  hist_met_diff->GetYaxis()->SetLabelSize(0.08);
  hist_met_diff->SetYTitle("#frac{DATA}{MC}");
  hist_met_diff->GetYaxis()->SetTitleSize(0.08);
  hist_met_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_met_diff->SetFillColorAlpha(45,0.35);
  hist_met_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// lead_pt ///
  TCanvas *c2 = new TCanvas("c2", "", 800, 800);
  TPad *c2_up = new TPad("c2_up", "", 0, 0.25, 1, 1);
  c2_up->SetTopMargin( 0.05 ); c2_up->SetBottomMargin( 0.02 ); c2_up->SetRightMargin( 0.02 ); c2_up->SetLeftMargin( 0.1 );
  TPad *c2_down = new TPad("c2_down", "", 0, 0, 1, 0.25);
  c2_down->SetTopMargin( 0.03 ); c2_down->SetBottomMargin( 0.25 ); c2_down->SetRightMargin( 0.02 ); c2_down->SetLeftMargin( 0.1 ); c2_down->SetGridx(); c2_down->SetGridy();
  c2_up->Draw();
  c2_down->Draw();
  c2_up->cd();
  ///bkg///
  bkgstack_lead_pt->Draw("hist");
  bkgstack_lead_pt->SetMaximum(y_max[1]);
  bkgstack_lead_pt->GetXaxis()->SetLabelSize(0);
  bkgstack_lead_pt->GetYaxis()->SetTitle("Events / "+TString::Itoa(lead_pt_dx,10)+" GeV");
  bkgstack_lead_pt->GetYaxis()->SetLabelSize(0.04);
  bkgstack_lead_pt->GetYaxis()->SetTitleSize(0.04);
  bkgstack_lead_pt->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_lead_pt[1]->SetMarkerStyle(3);
  hist_lead_pt[1]->SetMarkerSize(1);
  hist_lead_pt[1]->Draw("APsameE1");
  ///signal///
  hist_lead_pt[0]->SetLineColor(kRed);
  hist_lead_pt[0]->Draw("histsame");
  ///err///
  hist_lead_pt_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_lead_pt_err->SetFillStyle(3004);
  hist_lead_pt_err->SetFillColor(kBlue);
  hist_lead_pt_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c2_down->cd();
  TH1D *hist_lead_pt_diff = (TH1D*)hist_lead_pt[1]->Clone();
  hist_lead_pt_diff->Divide(hist_lead_pt_err);
  //hist_lead_pt_diff->Add(hist_lead_pt_err,-1);
  hist_lead_pt_diff->SetMaximum(2);
  hist_lead_pt_diff->SetMinimum(0);
  hist_lead_pt_diff->GetXaxis()->SetTitle("P_{T} of the leading lepton [GeV]");
  hist_lead_pt_diff->GetXaxis()->SetLabelSize(0.10);
  hist_lead_pt_diff->GetXaxis()->SetTitleSize(0.10);
  hist_lead_pt_diff->GetYaxis()->SetLabelSize(0.08);
  hist_lead_pt_diff->SetYTitle("#frac{DATA}{MC}");
  hist_lead_pt_diff->GetYaxis()->SetTitleSize(0.08);
  hist_lead_pt_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_lead_pt_diff->SetFillColorAlpha(45,0.35);
  hist_lead_pt_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// lead_eta ///
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  TPad *c3_up = new TPad("c3_up", "", 0, 0.25, 1, 1);
  c3_up->SetTopMargin( 0.05 ); c3_up->SetBottomMargin( 0.02 ); c3_up->SetRightMargin( 0.02 ); c3_up->SetLeftMargin( 0.1 );
  TPad *c3_down = new TPad("c3_down", "", 0, 0, 1, 0.25);
  c3_down->SetTopMargin( 0.03 ); c3_down->SetBottomMargin( 0.25 ); c3_down->SetRightMargin( 0.02 ); c3_down->SetLeftMargin( 0.1 ); c3_down->SetGridx(); c3_down->SetGridy();
  c3_up->Draw();
  c3_down->Draw();
  c3_up->cd();
  ///bkg///
  bkgstack_lead_eta->Draw("hist");
  bkgstack_lead_eta->SetMaximum(y_max[2]);
  bkgstack_lead_eta->GetXaxis()->SetLabelSize(0);
  bkgstack_lead_eta->GetYaxis()->SetTitle("Events / 0."+TString::Itoa(10*lead_eta_dx,10));
  bkgstack_lead_eta->GetYaxis()->SetLabelSize(0.04);
  bkgstack_lead_eta->GetYaxis()->SetTitleSize(0.04);
  bkgstack_lead_eta->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_lead_eta[1]->SetMarkerStyle(3);
  hist_lead_eta[1]->SetMarkerSize(1);
  hist_lead_eta[1]->Draw("APsameE1");
  ///signal///
  hist_lead_eta[0]->SetLineColor(kRed);
  hist_lead_eta[0]->Draw("histsame");
  ///err///
  hist_lead_eta_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_lead_eta_err->SetFillStyle(3004);
  hist_lead_eta_err->SetFillColor(kBlue);
  hist_lead_eta_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c3_down->cd();
  TH1D *hist_lead_eta_diff = (TH1D*)hist_lead_eta[1]->Clone();
  hist_lead_eta_diff->Divide(hist_lead_eta_err);
  //hist_lead_eta_diff->Add(hist_lead_eta_err,-1);
  hist_lead_eta_diff->SetMaximum(2);
  hist_lead_eta_diff->SetMinimum(0);
  hist_lead_eta_diff->GetXaxis()->SetTitle("#eta of the leading lepton");
  hist_lead_eta_diff->GetXaxis()->SetLabelSize(0.10);
  hist_lead_eta_diff->GetXaxis()->SetTitleSize(0.10);
  hist_lead_eta_diff->GetYaxis()->SetLabelSize(0.08);
  hist_lead_eta_diff->SetYTitle("#frac{DATA}{MC}");
  hist_lead_eta_diff->GetYaxis()->SetTitleSize(0.08);
  hist_lead_eta_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_lead_eta_diff->SetFillColorAlpha(45,0.35);
  hist_lead_eta_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// n_jets ///
  TCanvas *c4 = new TCanvas("c4", "", 800, 800);
  TPad *c4_up = new TPad("c4_up", "", 0, 0.25, 1, 1);
  c4_up->SetTopMargin( 0.05 ); c4_up->SetBottomMargin( 0.02 ); c4_up->SetRightMargin( 0.02 ); c4_up->SetLeftMargin( 0.1 );
  TPad *c4_down = new TPad("c4_down", "", 0, 0, 1, 0.25);
  c4_down->SetTopMargin( 0.03 ); c4_down->SetBottomMargin( 0.25 ); c4_down->SetRightMargin( 0.02 ); c4_down->SetLeftMargin( 0.1 ); c4_down->SetGridx(); c4_down->SetGridy();
  c4_up->Draw();
  c4_down->Draw();
  c4_up->cd();
  ///bkg///
  bkgstack_n_jets->Draw("hist");
  bkgstack_n_jets->SetMaximum(y_max[3]);
  bkgstack_n_jets->GetXaxis()->SetLabelSize(0);
  bkgstack_n_jets->GetYaxis()->SetTitle("Events");
  bkgstack_n_jets->GetYaxis()->SetLabelSize(0.04);
  bkgstack_n_jets->GetYaxis()->SetTitleSize(0.04);
  bkgstack_n_jets->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_n_jets[1]->SetMarkerStyle(3);
  hist_n_jets[1]->SetMarkerSize(1);
  hist_n_jets[1]->Draw("APsameE1");
  ///signal///
  hist_n_jets[0]->SetLineColor(kRed);
  hist_n_jets[0]->Draw("histsame");
  ///err///
  hist_n_jets_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_n_jets_err->SetFillStyle(3004);
  hist_n_jets_err->SetFillColor(kBlue);
  hist_n_jets_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c4_down->cd();
  TH1D *hist_n_jets_diff = (TH1D*)hist_n_jets[1]->Clone();
  hist_n_jets_diff->Divide(hist_n_jets_err);
  //hist_n_jets_diff->Add(hist_n_jets_err,-1);
  hist_n_jets_diff->SetMaximum(2);
  hist_n_jets_diff->SetMinimum(0);
  hist_n_jets_diff->GetXaxis()->SetTitle("# of jets");
  hist_n_jets_diff->GetXaxis()->SetLabelSize(0.10);
  hist_n_jets_diff->GetXaxis()->SetTitleSize(0.10);
  hist_n_jets_diff->GetYaxis()->SetLabelSize(0.08);
  hist_n_jets_diff->SetYTitle("#frac{DATA}{MC}");
  hist_n_jets_diff->GetYaxis()->SetTitleSize(0.08);
  hist_n_jets_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_n_jets_diff->SetFillColorAlpha(45,0.35);
  hist_n_jets_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// n_bjets ///
  TCanvas *c5 = new TCanvas("c5", "", 800, 800);
  TPad *c5_up = new TPad("c5_up", "", 0, 0.25, 1, 1);
  c5_up->SetTopMargin( 0.05 ); c5_up->SetBottomMargin( 0.02 ); c5_up->SetRightMargin( 0.02 ); c5_up->SetLeftMargin( 0.1 );
  TPad *c5_down = new TPad("c5_down", "", 0, 0, 1, 0.25);
  c5_down->SetTopMargin( 0.03 ); c5_down->SetBottomMargin( 0.25 ); c5_down->SetRightMargin( 0.02 ); c5_down->SetLeftMargin( 0.1 ); c5_down->SetGridx(); c5_down->SetGridy();
  c5_up->Draw();
  c5_down->Draw();
  c5_up->cd();
  c5_up->SetLogy();
  ///bkg///
  bkgstack_n_bjets->Draw("hist");
  bkgstack_n_bjets->SetMaximum(y_max[4]);
  bkgstack_n_bjets->GetXaxis()->SetLabelSize(0);
  bkgstack_n_bjets->GetYaxis()->SetTitle("Events");
  bkgstack_n_bjets->GetYaxis()->SetLabelSize(0.04);
  bkgstack_n_bjets->GetYaxis()->SetTitleSize(0.04);
  bkgstack_n_bjets->GetYaxis()->SetTitleOffset(1.2);
  ///data///
  hist_n_bjets[1]->SetMarkerStyle(3);
  hist_n_bjets[1]->SetMarkerSize(1);
  hist_n_bjets[1]->Draw("APsameE1");
  ///signal///
  hist_n_bjets[0]->SetLineColor(kRed);
  hist_n_bjets[0]->Draw("histsame");
  ///err///
  hist_n_bjets_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_n_bjets_err->SetFillStyle(3004);
  hist_n_bjets_err->SetFillColor(kBlue);
  hist_n_bjets_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c5_down->cd();
  TH1D *hist_n_bjets_diff = (TH1D*)hist_n_bjets[1]->Clone();
  hist_n_bjets_diff->Divide(hist_n_bjets_err);
  //hist_n_bjets_diff->Add(hist_n_bjets_err,-1);
  hist_n_bjets_diff->SetMaximum(2);
  hist_n_bjets_diff->SetMinimum(0);
  hist_n_bjets_diff->GetXaxis()->SetTitle("# of b-tagged jets");
  hist_n_bjets_diff->GetXaxis()->SetLabelSize(0.10);
  hist_n_bjets_diff->GetXaxis()->SetTitleSize(0.10);
  hist_n_bjets_diff->GetYaxis()->SetLabelSize(0.08);
  hist_n_bjets_diff->SetYTitle("#frac{DATA}{MC}");
  hist_n_bjets_diff->GetYaxis()->SetTitleSize(0.08);
  hist_n_bjets_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_n_bjets_diff->SetFillColorAlpha(45,0.35);
  hist_n_bjets_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  if(isSave){
    c1->SaveAs(savepath+TString::Itoa(sample, 10)+"_met.root");
    c2->SaveAs(savepath+TString::Itoa(sample, 10)+"_lead_pt.root");
    c3->SaveAs(savepath+TString::Itoa(sample, 10)+"_lead_eta.root");
    c4->SaveAs(savepath+TString::Itoa(sample, 10)+"_n_jets.root");
    c5->SaveAs(savepath+TString::Itoa(sample, 10)+"_n_bjets.root");
  }
  if(isClose){
    c1->Close();
    c2->Close();
    c3->Close();
    c4->Close();
    c5->Close();
  }
  cout
  << endl
  << "Logfile ------> " << logger_class.GetLogfilePath() << endl;
  logger_class.write_onebin(hist_sig, hist_dat, hist_bkg);
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
    
    MET_temp = MET_recon + lambda*MET_width;
    nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
    
    if(MET_temp >= 0){ // MET should be positive
      Pz_temp = solveqdeq(&(fitresult->isSamePz), 80.4, l1l2l3, MET_temp, METphi_recon, pm);
      PutNuPz(&nu, Pz_temp);
      W_real = l1l2l3 + nu;
      x_temp = chisquare(W_real.M(), lambda);
      if(W_real.M()==80.4) break; // lambda = 0 solvable
      if(x_temp < x){
        x = x_temp;
        fitresult->chi2 = x;
        fitresult->MET_fit = MET_temp;
        fitresult->lambda_fit = -range_lambda+dlambda*i;
        fitresult->Pz_fit = Pz_temp;
      }
    }
  }
}

int findleadlep(TLorentzVector *lep){
  double MaxPt = TMath::Max( lep[0].Pt() , TMath::Max( lep[1].Pt() , lep[2].Pt() ) );
  if(lep[0].Pt() == MaxPt) return 0;
  else if(lep[1].Pt() == MaxPt) return 1;
  else return 2;
}
