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

#define W_width 2.1
#define MET_width 7.0
#define N_chain 25
#define N_MCsector 7

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

void chi2_lambda_stack_mumumu(){
  
  TH1::SetDefaultSumw2(true);
  
  TStopwatch totaltime;
  totaltime.Start();
  int sample, isPU, ischi2, additional_cut, isSave=0, isClose=0, y_max[5];
  double HN_x_min, HN_x_max, HN_dx, W_real_x_min, W_real_x_max, W_real_dx, dR_x_min, dR_x_max, dR_dx, gamma_star_x_min, gamma_star_x_max, gamma_star_dx, z_candidate_x_min, z_candidate_x_max, z_candidate_dx;
  TString savepath = "/home/jskim/Documents/CMS/trilepton/mumumu/plots/";
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
  << "  22: delta R > 0.5 & W_reco < 100 GeV & HN Mass cut([30,80] for 40/50, [50,80] for 60)" << endl
  << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
  << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda = 0" << endl
  << "===============================================" << endl
  << "Cut : ";
  cin >> additional_cut;
  cout << endl;
  
  if(additional_cut==0){
    HN_x_min=0; HN_x_max=500; HN_dx=10; y_max[0]= 500;
    //W_real_x_min=80; W_real_x_max=81; W_real_dx=0.01; y_max[1] = 10;
		W_real_x_min=70; W_real_x_max=500; W_real_dx=10; y_max[1] = 600;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1; y_max[2] = 600;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1; y_max[3] = 600;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1; y_max[4] = 600;
    savepath += "before_dR/";
    ischi2 = 0;
  }
  if(additional_cut==1){
    HN_x_min=0; HN_x_max=500; HN_dx=10; y_max[0] = 350;
    W_real_x_min=70; W_real_x_max=500; W_real_dx=10; y_max[1] = 600;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1; y_max[2] = 200;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1; y_max[3] = 200;
    z_candidate_x_min=0, z_candidate_x_max=150, z_candidate_dx=10;y_max[4] = 2000;
    savepath += "dR/";
    ischi2 = 0;
  }
  if(additional_cut==2){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 100;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 100;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 60;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 60;
    savepath += "dR_W/";
    ischi2 = 0;
  }
  if(additional_cut==22){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 100;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 100;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 60;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 60;
    savepath += "dR_W_HN/";
    ischi2 = 0;
  }
  if(additional_cut==3){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 70;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 80;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 40;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 40;
    savepath += "dR_chi2/";
    ischi2 = 1;
  }
  if(additional_cut==4){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 70;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 80;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 40;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 40;
    savepath += "dR_chi2_lambda/";
    ischi2 = 1;
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
  
  TH1D *hist_HN[N_chain];
  TH1D *hist_W_real[N_chain];
  TH1D *hist_dR[N_chain];
  TH1D *hist_gamma_star[N_chain];
  TH1D *hist_z_candidate[N_chain];
  
  TH1D *hist_HN_err = new TH1D("", "", (HN_x_max-HN_x_min)/HN_dx, HN_x_min, HN_x_max); // 0, 100, 10
  TH1D *hist_W_real_err = new TH1D("", "", (W_real_x_max-W_real_x_min)/W_real_dx, W_real_x_min, W_real_x_max); // 70, 100, 30
  TH1D *hist_dR_err = new TH1D("", "", (dR_x_max-dR_x_min)/dR_dx, dR_x_min, dR_x_max);
  TH1D *hist_gamma_star_err = new TH1D("", "", (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx, gamma_star_x_min, gamma_star_x_max);
  TH1D *hist_z_candidate_err = new TH1D("", "", (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx, z_candidate_x_min, z_candidate_x_max);
  
  THStack *bkgstack_HN = new THStack("", "");
  THStack *bkgstack_W_real = new THStack("", "");
  THStack *bkgstack_dR = new THStack("", "");
  THStack *bkgstack_gamma_star = new THStack("", "");
  THStack *bkgstack_z_candidate = new THStack("", "");
  
  TH1D *hist_sig = new TH1D("N_sig", "N_sig", 1, 0, 1);
  TH1D *hist_dat = new TH1D("N_dat", "N_dat", 1, 0, 1);
  TH1D *hist_bkg = new TH1D("N_bkg", "N_bkg", 1, 0, 1);
  TH1D *hist_sec[N_MCsector];
  
  TLorentzVector lep[3], W_real, HN, nu, gen_nu, gamma_star, z_candidate, selection_nu[2];

  int FillOrder[N_chain] = { 0, 1, // signal, data
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, // others
    22, 23, // Higgs
    24, // Wtollln
    5, 6, 7, 8, 9, 10, 11, // VV
    2, 3, 4 }, sec_num = 0; // DY
  TString MCsector[N_MCsector] = {"signal", "data", "others", "Higgs", "Wtollln", "VV", "DY"};
  
  for(int i=0;i<N_MCsector;i++){
    hist_sec[i] = new TH1D("N_"+MCsector[i], "N_"+MCsector[i], 1, 0, 1);
  }
  
  for(int l=0; l<N_chain; l++){
    int k = FillOrder[l];
    N_entry = 0;
    if(fillcolors[FillOrder[l-1]] != fillcolors[FillOrder[l]] ){ // sector changed
      logger_class.save_n_entry_sec(hist_sec[sec_num]->GetBinContent(1), hist_sec[sec_num]->GetBinError(1));
      sec_num++;
    }
    if(l==N_chain-1){ // last sector
      logger_class.save_n_entry_sec(hist_sec[sec_num]->GetBinContent(1), hist_sec[sec_num]->GetBinError(1));
    }
    PUweight_sum = 0;
 
    mu[k] = new Muon(chain_mu[k], 3);
    evt[k] = new Event(chain_event[k]);
    jet[k] = new Jet(chain_jet[k], 4);
    
    hist_HN[k] = new TH1D("", "", (HN_x_max-HN_x_min)/HN_dx, HN_x_min, HN_x_max);
    hist_W_real[k] = new TH1D("", "", (W_real_x_max-W_real_x_min)/W_real_dx, W_real_x_min, W_real_x_max);
    hist_dR[k] = new TH1D("", "", (dR_x_max-dR_x_min)/dR_dx, dR_x_min, dR_x_max);
    hist_gamma_star[k] = new TH1D("", "", (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx, gamma_star_x_min, gamma_star_x_max);
    hist_z_candidate[k] = new TH1D("", "", (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx, z_candidate_x_min, z_candidate_x_max);

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
      selection = mu[k]->isHNLooseMuon(0, 15) && mu[k]->isHNLooseMuon(1, 10) && mu[k]->isHNLooseMuon(2, 10);
      //selection = mu[k]->isHNTightMuon(0, 15) && mu[k]->isHNTightMuon(1, 10) && mu[k]->isHNTightMuon(2, 10);
      //Cut : dR > 0.5//
      if(additional_cut == 1 || additional_cut == 2  || additional_cut == 22 || additional_cut == 3 || additional_cut == 4){
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
        if(additional_cut == 2 || additional_cut == 22){
          selection *= W_real.M() < 100.0;
        }
        //Cut : HN(Reco) cut //
        if(additional_cut == 22) {
          if(sample == 40 || sample == 50){
            selection *= HN.M() > 30 && HN.M() < 80;
          }
          if(sample == 60){
            selection *= HN.M() > 50 && HN.M() < 80;
          }
        }
        //Cut : chi2 < 0.1//
        if(additional_cut == 3 || additional_cut == 3){
          selection *= fitresult[solution_selection].chi2 < 0.1;
        }
        //Cut : lambda = 0//
        if(additional_cut == 4){
          selection *= fitresult[solution_selection].lambda_fit == 0;
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
          
          hist_sec[sec_num]->Fill(0., weight);
          
          // histogram filling //
          hist_HN[k]->Fill(HN.M(), weight);
          hist_W_real[k]->Fill(W_real.M(), weight);
          hist_dR[k]->Fill( TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ) , weight);
          hist_gamma_star[k]->Fill(gamma_star.M(), weight);
          hist_z_candidate[k]->Fill(z_candidate.M(), weight);
          
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
    
    hist_HN[k]->SetLineColor(fillcolors[k]);
    hist_W_real[k]->SetLineColor(fillcolors[k]);
    hist_dR[k]->SetLineColor(fillcolors[k]);
    hist_gamma_star[k]->SetLineColor(fillcolors[k]);
    hist_z_candidate[k]->SetLineColor(fillcolors[k]);
    hist_HN[k]->SetFillColor(fillcolors[k]);
    hist_W_real[k]->SetFillColor(fillcolors[k]);
    hist_dR[k]->SetFillColor(fillcolors[k]);
    hist_gamma_star[k]->SetFillColor(fillcolors[k]);
    hist_z_candidate[k]->SetFillColor(fillcolors[k]);
    
    /// make TH1 for error & make THSatck
    if(k>=2){
      bkgstack_HN->Add(hist_HN[k]);
      bkgstack_W_real->Add(hist_W_real[k]);
      bkgstack_dR->Add(hist_dR[k]);
      bkgstack_gamma_star->Add(hist_gamma_star[k]);
      bkgstack_z_candidate->Add(hist_z_candidate[k]);
      hist_HN_err->Add(hist_HN[k]);
      hist_W_real_err->Add(hist_W_real[k]);
      hist_dR_err->Add(hist_dR[k]);
      hist_gamma_star_err->Add(hist_gamma_star[k]);
      hist_z_candidate_err->Add(hist_z_candidate[k]);
    }

  
  } // chain loop
  
  logger_class.divide_log();
  vector<TString> v_MCsector(MCsector, MCsector+N_MCsector);
  logger_class.write_n_entry_sec(v_MCsector);

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
  hist_HN[0]->Scale(k_factor*coupling_const);
  hist_W_real[0]->Scale(k_factor*coupling_const);
  hist_dR[0]->Scale(k_factor*coupling_const);
  hist_gamma_star[0]->Scale(k_factor*coupling_const);
  hist_z_candidate[0]->Scale(k_factor*coupling_const);
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  //legend//
  TLegend *lg1 = new TLegend(0.65, 0.55, 0.97, 0.90);
  lg1->SetFillStyle(0);
  lg1->SetBorderSize(0);
  lg1->AddEntry(hist_HN[1], "data", "p");
  lg1->AddEntry(hist_HN[2], "DY", "f");
  lg1->AddEntry(hist_HN[5], "VV", "f");
  lg1->AddEntry(hist_HN[24], "W#rightarrowlll#nu", "f");
  lg1->AddEntry(hist_HN[22], "Higgs", "f");
  lg1->AddEntry(hist_HN[14], "others", "f");
  lg1->AddEntry(hist_HN[0], "HN"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  //y=0//
  double x_0[2], y_0[2];
  x_0[0] = -1000;  y_0[0] = 0;
  x_0[1] = 1000;  y_0[1] = 0;
  TGraph *g0 = new TGraph(2, x_0, y_0);
  //y=1//
  double x_1[2], y_1[2];
  x_1[0] = 1000;  y_1[0] = 1;
  x_1[1] = -1000;  y_1[1] = 1;
  TGraph *g1 = new TGraph(2, x_1, y_1);
  
  /// HN ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  TPad *c1_up = new TPad("c1_up", "", 0, 0.25, 1, 1);
  c1_up->SetTopMargin( 0.05 ); c1_up->SetBottomMargin( 0.02 ); c1_up->SetRightMargin( 0.02 ); c1_up->SetLeftMargin( 0.1 );
  TPad *c1_down = new TPad("c1_down", "", 0, 0, 1, 0.25);
  c1_down->SetTopMargin( 0.03 ); c1_down->SetBottomMargin( 0.25 ); c1_down->SetRightMargin( 0.02 ); c1_down->SetLeftMargin( 0.1 ); c1_down->SetGridx(); c1_down->SetGridy();
  c1_up->Draw();
  c1_down->Draw();
  c1_up->cd();
  ///bkg///
  bkgstack_HN->Draw("hist");
  bkgstack_HN->SetMaximum(y_max[0]);
  bkgstack_HN->GetXaxis()->SetLabelSize(0);
  bkgstack_HN->GetYaxis()->SetTitle("Events / "+TString::Itoa(HN_dx,10)+" GeV");
  bkgstack_HN->GetYaxis()->SetLabelSize(0.05);
  bkgstack_HN->GetYaxis()->SetTitleSize(0.05);
  bkgstack_HN->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_HN[1]->SetMarkerStyle(3);
  hist_HN[1]->SetMarkerSize(1);
  hist_HN[1]->Draw("APsameE1");
  ///signal///
  hist_HN[0]->SetLineColor(kRed);
  hist_HN[0]->Draw("histsame");
  ///err///
  hist_HN_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_HN_err->SetFillStyle(3004);
  hist_HN_err->SetFillColor(kBlue);
  hist_HN_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c1_down->cd();
  TH1D *hist_HN_diff = (TH1D*)hist_HN[1]->Clone();
  hist_HN_diff->Divide(hist_HN_err);
  //hist_HN_diff->Add(hist_HN_err,-1);
  hist_HN_diff->SetMaximum(2);
  hist_HN_diff->SetMinimum(0);
  hist_HN_diff->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  hist_HN_diff->GetXaxis()->SetLabelSize(0.10);
  hist_HN_diff->GetXaxis()->SetTitleSize(0.10);
  hist_HN_diff->GetYaxis()->SetLabelSize(0.08);
  hist_HN_diff->SetYTitle("#frac{DATA}{MC}");
  hist_HN_diff->GetYaxis()->SetTitleSize(0.08);
  hist_HN_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_HN_diff->SetFillColorAlpha(45,0.35);
  hist_HN_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");

  /// W_real ///
  TCanvas *c2 = new TCanvas("c2", "", 800, 800);
  TPad *c2_up = new TPad("c2_up", "", 0, 0.25, 1, 1);
  c2_up->SetTopMargin( 0.05 ); c2_up->SetBottomMargin( 0.02 ); c2_up->SetRightMargin( 0.02 ); c2_up->SetLeftMargin( 0.1 );
  TPad *c2_down = new TPad("c2_down", "", 0, 0, 1, 0.25);
  c2_down->SetTopMargin( 0.03 ); c2_down->SetBottomMargin( 0.25 ); c2_down->SetRightMargin( 0.02 ); c2_down->SetLeftMargin( 0.1 ); c2_down->SetGridx(); c2_down->SetGridy();
  c2_up->Draw();
  c2_down->Draw();
  c2_up->cd();
  ///bkg///
  bkgstack_W_real->Draw("hist");
  bkgstack_W_real->SetMaximum(y_max[1]);
  bkgstack_W_real->GetXaxis()->SetLabelSize(0);
  bkgstack_W_real->GetYaxis()->SetTitle("Events / "+TString::Itoa(W_real_dx,10)+" GeV");
  bkgstack_W_real->GetYaxis()->SetLabelSize(0.05);
  bkgstack_W_real->GetYaxis()->SetTitleSize(0.05);
  bkgstack_W_real->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_W_real[1]->SetMarkerStyle(3);
  hist_W_real[1]->SetMarkerSize(1);
  hist_W_real[1]->Draw("APsameE1");
  ///signal///
  hist_W_real[0]->SetLineColor(kRed);
  hist_W_real[0]->Draw("histsame");
  ///err///
  hist_W_real_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_W_real_err->SetFillStyle(3004);
  hist_W_real_err->SetFillColor(kBlue);
  hist_W_real_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c2_down->cd();
  TH1D *hist_W_real_diff = (TH1D*)hist_W_real[1]->Clone();
  hist_W_real_diff->Divide(hist_W_real_err);
  //hist_W_real_diff->Add(hist_W_real_err,-1);
  hist_W_real_diff->SetMaximum(2);
  hist_W_real_diff->SetMinimum(0);
  hist_W_real_diff->GetXaxis()->SetTitle("m(#mu_{1}#mu_{2}#mu_{3}#nu) [GeV]");
  hist_W_real_diff->GetXaxis()->SetLabelSize(0.10);
  hist_W_real_diff->GetXaxis()->SetTitleSize(0.10);
  hist_W_real_diff->GetYaxis()->SetLabelSize(0.08);
  hist_W_real_diff->SetYTitle("#frac{DATA}{MC}");
  hist_W_real_diff->GetYaxis()->SetTitleSize(0.08);
  hist_W_real_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_W_real_diff->SetFillColorAlpha(45,0.35);
  hist_W_real_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");

  /// dR ///
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  TPad *c3_up = new TPad("c3_up", "", 0, 0.25, 1, 1);
  c3_up->SetTopMargin( 0.05 ); c3_up->SetBottomMargin( 0.02 ); c3_up->SetRightMargin( 0.02 ); c3_up->SetLeftMargin( 0.1 );
  TPad *c3_down = new TPad("c3_down", "", 0, 0, 1, 0.25);
  c3_down->SetTopMargin( 0.03 ); c3_down->SetBottomMargin( 0.25 ); c3_down->SetRightMargin( 0.02 ); c3_down->SetLeftMargin( 0.1 ); c3_down->SetGridx(); c3_down->SetGridy();
  c3_up->Draw();
  c3_down->Draw();
  c3_up->cd();
  ///bkg///
  bkgstack_dR->Draw("hist");
  bkgstack_dR->SetMaximum(y_max[2]);
  bkgstack_dR->GetXaxis()->SetLabelSize(0);
  bkgstack_dR->GetYaxis()->SetTitle("Events / 0."+TString::Itoa(10*dR_dx,10));
  bkgstack_dR->GetYaxis()->SetLabelSize(0.05);
  bkgstack_dR->GetYaxis()->SetTitleSize(0.05);
  bkgstack_dR->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_dR[1]->SetMarkerStyle(3);
  hist_dR[1]->SetMarkerSize(1);
  hist_dR[1]->Draw("APsameE1");
  ///signal///
  hist_dR[0]->SetLineColor(kRed);
  hist_dR[0]->Draw("histsame");
  ///err///
  hist_dR_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_dR_err->SetFillStyle(3004);
  hist_dR_err->SetFillColor(kBlue);
  hist_dR_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c3_down->cd();
  TH1D *hist_dR_diff = (TH1D*)hist_dR[1]->Clone();
  hist_dR_diff->Divide(hist_dR_err);
  //hist_dR_diff->Add(hist_dR_err,-1);
  hist_dR_diff->SetMaximum(2);
  hist_dR_diff->SetMinimum(0);
  hist_dR_diff->GetXaxis()->SetTitle("#DeltaR(OS)_{Min}");
  hist_dR_diff->GetXaxis()->SetLabelSize(0.10);
  hist_dR_diff->GetXaxis()->SetTitleSize(0.10);
  hist_dR_diff->GetYaxis()->SetLabelSize(0.08);
  hist_dR_diff->SetYTitle("#frac{DATA}{MC}");
  hist_dR_diff->GetYaxis()->SetTitleSize(0.08);
  hist_dR_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_dR_diff->SetFillColorAlpha(45,0.35);
  hist_dR_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
 
  /// gamma_star ///
  TCanvas *c4 = new TCanvas("c4", "", 800, 800);
  TPad *c4_up = new TPad("c4_up", "", 0, 0.25, 1, 1);
  c4_up->SetTopMargin( 0.05 ); c4_up->SetBottomMargin( 0.02 ); c4_up->SetRightMargin( 0.02 ); c4_up->SetLeftMargin( 0.1 );
  TPad *c4_down = new TPad("c4_down", "", 0, 0, 1, 0.25);
  c4_down->SetTopMargin( 0.03 ); c4_down->SetBottomMargin( 0.25 ); c4_down->SetRightMargin( 0.02 ); c4_down->SetLeftMargin( 0.1 ); c4_down->SetGridx(); c4_down->SetGridy();
  c4_up->Draw();
  c4_down->Draw();
  c4_up->cd();
  ///bkg///
  bkgstack_gamma_star->Draw("hist");
  bkgstack_gamma_star->SetMaximum(y_max[3]);
  bkgstack_gamma_star->GetXaxis()->SetLabelSize(0);
  bkgstack_gamma_star->GetYaxis()->SetTitle("Events / "+TString::Itoa(gamma_star_dx,10)+" GeV");
  bkgstack_gamma_star->GetYaxis()->SetLabelSize(0.05);
  bkgstack_gamma_star->GetYaxis()->SetTitleSize(0.05);
  bkgstack_gamma_star->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_gamma_star[1]->SetMarkerStyle(3);
  hist_gamma_star[1]->SetMarkerSize(1);
  hist_gamma_star[1]->Draw("APsameE1");
  ///signal///
  hist_gamma_star[0]->SetLineColor(kRed);
  hist_gamma_star[0]->Draw("histsame");
  ///err///
  hist_gamma_star_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_gamma_star_err->SetFillStyle(3004);
  hist_gamma_star_err->SetFillColor(kBlue);
  hist_gamma_star_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c4_down->cd();
  TH1D *hist_gamma_star_diff = (TH1D*)hist_gamma_star[1]->Clone();
  hist_gamma_star_diff->Divide(hist_gamma_star_err);
  //hist_gamma_star_diff->Add(hist_gamma_star_err,-1);
  hist_gamma_star_diff->SetMaximum(2);
  hist_gamma_star_diff->SetMinimum(0);
  hist_gamma_star_diff->GetXaxis()->SetTitle("m(#mu+#mu-) for Minimum #DeltaR [GeV]");
  hist_gamma_star_diff->GetXaxis()->SetLabelSize(0.10);
  hist_gamma_star_diff->GetXaxis()->SetTitleSize(0.10);
  hist_gamma_star_diff->GetYaxis()->SetLabelSize(0.08);
  hist_gamma_star_diff->SetYTitle("#frac{DATA}{MC}");
  hist_gamma_star_diff->GetYaxis()->SetTitleSize(0.08);
  hist_gamma_star_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_gamma_star_diff->SetFillColorAlpha(45,0.35);
  hist_gamma_star_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");

  /// z_candidate ///
  TCanvas *c5 = new TCanvas("c5", "", 800, 800);
  TPad *c5_up = new TPad("c5_up", "", 0, 0.25, 1, 1);
  c5_up->SetTopMargin( 0.05 ); c5_up->SetBottomMargin( 0.02 ); c5_up->SetRightMargin( 0.02 ); c5_up->SetLeftMargin( 0.1 );
  TPad *c5_down = new TPad("c5_down", "", 0, 0, 1, 0.25);
  c5_down->SetTopMargin( 0.03 ); c5_down->SetBottomMargin( 0.25 ); c5_down->SetRightMargin( 0.02 ); c5_down->SetLeftMargin( 0.1 ); c5_down->SetGridx(); c5_down->SetGridy();
  c5_up->Draw();
  c5_down->Draw();
  c5_up->cd();
  ///bkg///
  bkgstack_z_candidate->Draw("hist");
  bkgstack_z_candidate->SetMaximum(y_max[4]);
  bkgstack_z_candidate->GetXaxis()->SetLabelSize(0);
  bkgstack_z_candidate->GetYaxis()->SetTitle("Events / 10 GeV");
  bkgstack_z_candidate->GetYaxis()->SetLabelSize(0.05);
  bkgstack_z_candidate->GetYaxis()->SetTitleSize(0.05);
  bkgstack_z_candidate->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_z_candidate[1]->SetMarkerStyle(3);
  hist_z_candidate[1]->SetMarkerSize(1);
  hist_z_candidate[1]->Draw("APsameE1");
  ///signal///
  hist_z_candidate[0]->SetLineColor(kRed);
  hist_z_candidate[0]->Draw("histsame");
  ///err///
  hist_z_candidate_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_z_candidate_err->SetFillStyle(3004);
  hist_z_candidate_err->SetFillColor(kBlue);
  hist_z_candidate_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c5_down->cd();
  TH1D *hist_z_candidate_diff = (TH1D*)hist_z_candidate[1]->Clone();
  hist_z_candidate_diff->Divide(hist_z_candidate_err);
  //hist_z_candidate_diff->Add(hist_z_candidate_err,-1);
  hist_z_candidate_diff->SetMaximum(2);
  hist_z_candidate_diff->SetMinimum(0);
  hist_z_candidate_diff->GetXaxis()->SetTitle("m(#mu+#mu-) for Z-candidate [GeV]");
  hist_z_candidate_diff->GetXaxis()->SetLabelSize(0.10);
  hist_z_candidate_diff->GetXaxis()->SetTitleSize(0.10);
  hist_z_candidate_diff->GetYaxis()->SetLabelSize(0.08);
  hist_z_candidate_diff->SetYTitle("#frac{DATA}{MC}");
  hist_z_candidate_diff->GetYaxis()->SetTitleSize(0.08);
  hist_z_candidate_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_z_candidate_diff->SetFillColorAlpha(45,0.35);
  hist_z_candidate_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");

  if(isSave){
    c1->SaveAs(savepath+TString::Itoa(sample, 10)+"_HN.pdf");
    c2->SaveAs(savepath+TString::Itoa(sample, 10)+"_W.pdf");
    c3->SaveAs(savepath+TString::Itoa(sample, 10)+"_dR.pdf");
    c4->SaveAs(savepath+TString::Itoa(sample, 10)+"_gamma_star.pdf");
    c5->SaveAs(savepath+TString::Itoa(sample, 10)+"_z_candidate.pdf");
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

void chi2_lambda_stack_mumumu(int sample, int isPU, int additional_cut, int muonidsel_int){ // muonidsel = 0 = loose, muonidsel = 1 = tight
  
  TH1::SetDefaultSumw2(true);
  
  TStopwatch totaltime;
  totaltime.Start();
  //int sample, isPU, ischi2, additional_cut;
  int isSave=1, isClose=1, y_max[5], ischi2;
  double HN_x_min, HN_x_max, HN_dx, W_real_x_min, W_real_x_max, W_real_dx, dR_x_min, dR_x_max, dR_dx, gamma_star_x_min, gamma_star_x_max, gamma_star_dx, z_candidate_x_min, z_candidate_x_max, z_candidate_dx;
  TString savepath = "/home/jskim/Documents/CMS/trilepton/mumumu/plots/", muonidsel;
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
   << "  22: delta R > 0.5 & W_reco < 100 GeV & HN Mass cut([30,80] for 40/50, [50,80] for 60)" << endl
   << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
   << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda = 0" << endl
   << "===============================================" << endl
   << "Cut : ";
   cin >> additional_cut;
   cout << endl;
   */
  
  if(muonidsel_int == 0)  muonidsel = "loose";
  if(muonidsel_int == 1)  muonidsel = "tight";
  
  if(additional_cut==0){
    HN_x_min=0; HN_x_max=500; HN_dx=10; y_max[0]= 500;
    W_real_x_min=70; W_real_x_max=500; W_real_dx=10; y_max[1] = 600;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1; y_max[2] = 600;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1; y_max[3] = 600;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1; y_max[4] = 600;
    if(muonidsel == "loose")  savepath += "before_dR/loose_";
    if(muonidsel == "tight")  savepath += "before_dR/tight_";
    ischi2 = 0;
  }
  if(additional_cut==1){
    HN_x_min=0; HN_x_max=500; HN_dx=10; y_max[0] = 350;
    W_real_x_min=70; W_real_x_max=500; W_real_dx=10; y_max[1] = 600;
    dR_x_min=0, dR_x_max=5, dR_dx=0.1; y_max[2] = 200;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=1; y_max[3] = 200;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=1;y_max[4] = 600;
    if(muonidsel == "loose")  savepath += "dR/loose_";
    if(muonidsel == "tight")  savepath += "dR/tight_";
    ischi2 = 0;
  }
  if(additional_cut==2){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 100;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 100;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 60;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 60;
    if(muonidsel == "loose")  savepath += "dR_W/loose_";
    if(muonidsel == "tight")  savepath += "dR_W/tight_";
    ischi2 = 0;
  }
  if(additional_cut==22){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 100;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 100;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 60;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 60;
    if(muonidsel == "loose")  savepath += "dR_W_HN/loose_";
    if(muonidsel == "tight")  savepath += "dR_W_HN/tight_";
    ischi2 = 0;
  }
  if(additional_cut==3){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 70;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 80;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 40;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 40;
    if(muonidsel == "loose")  savepath += "dR_chi2/loose_";
    if(muonidsel == "tight")  savepath += "dR_chi2/tight_";
    ischi2 = 1;
  }
  if(additional_cut==4){
    HN_x_min=0; HN_x_max=100; HN_dx=10; y_max[0] = 70;
    W_real_x_min=70; W_real_x_max=120; W_real_dx=5; y_max[1] = 100;
    dR_x_min=0, dR_x_max=5, dR_dx=0.5; y_max[2] = 80;
    gamma_star_x_min=0, gamma_star_x_max=120, gamma_star_dx=5; y_max[3] = 40;
    z_candidate_x_min=0, z_candidate_x_max=120, z_candidate_dx=5; y_max[4] = 40;
    if(muonidsel == "loose")  savepath += "dR_chi2_lambda/loose_";
    if(muonidsel == "tight")  savepath += "dR_chi2_lambda/tight_";
    ischi2 = 1;
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
  int N_event, solution_selection, smaller;
  struct Fit fitresult[2];
  Color_t fillcolors[N_chain] = {0, kBlack, // signal and data
    kAzure+8, kAzure+8, kAzure+8, // DY
    kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, kGreen, // VV
    kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, kRed, // others
    kViolet, kViolet, // higgs
    kYellow}; // Wtollln
  
  TH1D *hist_HN[N_chain];
  TH1D *hist_W_real[N_chain];
  TH1D *hist_dR[N_chain];
  TH1D *hist_gamma_star[N_chain];
  TH1D *hist_z_candidate[N_chain];
  
  TH1D *hist_HN_err = new TH1D("hist_HN_err", "", (HN_x_max-HN_x_min)/HN_dx, HN_x_min, HN_x_max);
  TH1D *hist_W_real_err = new TH1D("hist_W_real_err", "", (W_real_x_max-W_real_x_min)/W_real_dx, W_real_x_min, W_real_x_max);
  TH1D *hist_dR_err = new TH1D("hist_dR_err", "", (dR_x_max-dR_x_min)/dR_dx, dR_x_min, dR_x_max);
  TH1D *hist_gamma_star_err = new TH1D("hist_gamma_star_err", "", (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx, gamma_star_x_min, gamma_star_x_max);
  TH1D *hist_z_candidate_err = new TH1D("hist_z_candidate_err", "", (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx, z_candidate_x_min, z_candidate_x_max);
  
  THStack *bkgstack_HN = new THStack("bkgstack_HN", "");
  THStack *bkgstack_W_real = new THStack("bkgstack_W_real", "");
  THStack *bkgstack_dR = new THStack("bkgstack_dR", "");
  THStack *bkgstack_gamma_star = new THStack("bkgstack_gamma_star", "");
  THStack *bkgstack_z_candidate = new THStack("bkgstack_z_candidate", "");
  
  TH1D *hist_sig = new TH1D("N_sig", "N_sig", 1, 0, 1);
  TH1D *hist_dat = new TH1D("N_dat", "N_dat", 1, 0, 1);
  TH1D *hist_bkg = new TH1D("N_bkg", "N_bkg", 1, 0, 1);
  TH1D *hist_sec[N_MCsector];
  
  TLorentzVector lep[3], W_real, HN, nu, gen_nu, gamma_star, z_candidate, selection_nu[2];
  
  int FillOrder[N_chain] = { 0, 1, // signal, data
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, // others
    22, 23, // Higgs
    24, // Wtollln
    5, 6, 7, 8, 9, 10, 11, // VV
    2, 3, 4 }, sec_num = 0; // DY
  TString MCsector[N_MCsector] = {"signal", "data", "others", "Higgs", "Wtollln", "VV", "DY"};
  
  for(int i=0;i<N_MCsector;i++){
    hist_sec[i] = new TH1D("N_"+MCsector[i], "N_"+MCsector[i], 1, 0, 1);
  }
  
  for(int l=0; l<N_chain; l++){
    int k = FillOrder[l];
    N_event = 0;
    if(fillcolors[FillOrder[l-1]] != fillcolors[FillOrder[l]] ){ // sector changed
      logger_class.save_n_entry_sec(hist_sec[sec_num]->GetBinContent(1), hist_sec[sec_num]->GetBinError(1));
      sec_num++;
    }
    if(l==N_chain-1){ // last sector
      logger_class.save_n_entry_sec(hist_sec[sec_num]->GetBinContent(1), hist_sec[sec_num]->GetBinError(1));
    }
    PUweight_sum = 0;
    
    mu[k] = new Muon(chain_mu[k], 3);
    evt[k] = new Event(chain_event[k]);
    jet[k] = new Jet(chain_jet[k], 4);
    
    hist_HN[k] = new TH1D("hist_HN"+TString::Itoa(k,10), "", (HN_x_max-HN_x_min)/HN_dx, HN_x_min, HN_x_max);
    hist_W_real[k] = new TH1D("hist_W_real"+TString::Itoa(k,10), "", (W_real_x_max-W_real_x_min)/W_real_dx, W_real_x_min, W_real_x_max);
    hist_dR[k] = new TH1D("hist_dR"+TString::Itoa(k,10), "", (dR_x_max-dR_x_min)/dR_dx, dR_x_min, dR_x_max);
    hist_gamma_star[k] = new TH1D("hist_gamma_star"+TString::Itoa(k,10), "", (gamma_star_x_max-gamma_star_x_min)/gamma_star_dx, gamma_star_x_min, gamma_star_x_max);
    hist_z_candidate[k] = new TH1D("hist_z_candidate"+TString::Itoa(k,10), "", (z_candidate_x_max-z_candidate_x_min)/z_candidate_dx, z_candidate_x_min, z_candidate_x_max);
    
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
      if(additional_cut == 1 || additional_cut == 2  || additional_cut == 22 || additional_cut == 3 || additional_cut == 4){
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
        if(additional_cut == 2 || additional_cut == 22){
          selection *= W_real.M() < 100.0;
        }
        //Cut : HN(Reco) cut //
        if(additional_cut == 22) {
          if(sample == 40 || sample == 50){
            selection *= HN.M() > 30 && HN.M() < 80;
          }
          if(sample == 60){
            selection *= HN.M() > 50 && HN.M() < 80;
          }
        }
        //Cut : chi2 < 0.1//
        if(additional_cut == 3 || additional_cut == 3){
          selection *= fitresult[solution_selection].chi2 < 0.1;
        }
        //Cut : lambda = 0//
        if(additional_cut == 4){
          selection *= fitresult[solution_selection].lambda_fit == 0;
        }
        
        if(selection){
          N_event++; // pass-cut calculation
          
          // gamma_star, z_candidate 4-vec reconstruction //
          if( lep[0].DeltaR(lep[1]) < lep[2].DeltaR(lep[1]) ) gamma_star = lep[1] + lep[0];
          else                                                gamma_star = lep[1] + lep[2];
          if( fabs( (lep[1] + lep[0]).M()-91.19 ) >= fabs( (lep[1] + lep[2]).M()-91.19 ) ) z_candidate = lep[1] + lep[2];
          else                                                                             z_candidate = lep[1] + lep[0];
          
          // weight, PU total effect //
          if(!isPU) weight = evt[k]->weight;
          else      weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
          PUweight_sum+=evt[k]->PU_reweight;
          
          hist_sec[sec_num]->Fill(0., weight);
          
          // histogram filling //
          hist_HN[k]->Fill(HN.M(), weight);
          hist_W_real[k]->Fill(W_real.M(), weight);
          hist_dR[k]->Fill( TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ) , weight);
          hist_gamma_star[k]->Fill(gamma_star.M(), weight);
          hist_z_candidate[k]->Fill(z_candidate.M(), weight);
          
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
    << "event = " << N_event << endl
    << "PUweight = " << PUweight_sum/(double)N_event << endl;
    
    logger_class.write_n_events(chainlist[k], evt[k]->weight, N_event, PUweight_sum/(double)N_event);
    
    hist_HN[k]->SetLineColor(fillcolors[k]);
    hist_W_real[k]->SetLineColor(fillcolors[k]);
    hist_dR[k]->SetLineColor(fillcolors[k]);
    hist_gamma_star[k]->SetLineColor(fillcolors[k]);
    hist_z_candidate[k]->SetLineColor(fillcolors[k]);
    hist_HN[k]->SetFillColor(fillcolors[k]);
    hist_W_real[k]->SetFillColor(fillcolors[k]);
    hist_dR[k]->SetFillColor(fillcolors[k]);
    hist_gamma_star[k]->SetFillColor(fillcolors[k]);
    hist_z_candidate[k]->SetFillColor(fillcolors[k]);
    
    /// make TH1 for error
    if(k>=2){
      bkgstack_HN->Add(hist_HN[k]);
      bkgstack_W_real->Add(hist_W_real[k]);
      bkgstack_dR->Add(hist_dR[k]);
      bkgstack_gamma_star->Add(hist_gamma_star[k]);
      bkgstack_z_candidate->Add(hist_z_candidate[k]);
      hist_HN_err->Add(hist_HN[k]);
      hist_W_real_err->Add(hist_W_real[k]);
      hist_dR_err->Add(hist_dR[k]);
      hist_gamma_star_err->Add(hist_gamma_star[k]);
      hist_z_candidate_err->Add(hist_z_candidate[k]);
    }
    
    
  } // chain loop
  
  logger_class.divide_log();
  vector<TString> v_MCsector(MCsector, MCsector+N_MCsector);
  logger_class.write_n_entry_sec(v_MCsector);
  
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
  hist_HN[0]->Scale(k_factor*coupling_const);
  hist_W_real[0]->Scale(k_factor*coupling_const);
  hist_dR[0]->Scale(k_factor*coupling_const);
  hist_gamma_star[0]->Scale(k_factor*coupling_const);
  hist_z_candidate[0]->Scale(k_factor*coupling_const);
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  //legend//
  TLegend *lg1 = new TLegend(0.65, 0.55, 0.97, 0.90);
  lg1->SetFillStyle(0);
  lg1->SetBorderSize(0);
  lg1->AddEntry(hist_HN[1], "data", "p");
  lg1->AddEntry(hist_HN[2], "DY", "f");
  lg1->AddEntry(hist_HN[5], "VV", "f");
  lg1->AddEntry(hist_HN[24], "W#rightarrowlll#nu", "f");
  lg1->AddEntry(hist_HN[22], "Higgs", "f");
  lg1->AddEntry(hist_HN[14], "others", "f");
  lg1->AddEntry(hist_HN[0], "HN"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  //y=0//
  double x_0[2], y_0[2];
  x_0[0] = -1000;  y_0[0] = 0;
  x_0[1] = 1000;  y_0[1] = 0;
  TGraph *g0 = new TGraph(2, x_0, y_0);
  //y=1//
  double x_1[2], y_1[2];
  x_1[0] = 1000;  y_1[0] = 1;
  x_1[1] = -1000;  y_1[1] = 1;
  TGraph *g1 = new TGraph(2, x_1, y_1);
  
  /// HN ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 800);
  TPad *c1_up = new TPad("c1_up", "", 0, 0.25, 1, 1);
  c1_up->SetTopMargin( 0.05 ); c1_up->SetBottomMargin( 0.02 ); c1_up->SetRightMargin( 0.02 ); c1_up->SetLeftMargin( 0.1 );
  TPad *c1_down = new TPad("c1_down", "", 0, 0, 1, 0.25);
  c1_down->SetTopMargin( 0.03 ); c1_down->SetBottomMargin( 0.25 ); c1_down->SetRightMargin( 0.02 ); c1_down->SetLeftMargin( 0.1 ); c1_down->SetGridx(); c1_down->SetGridy();
  c1_up->Draw();
  c1_down->Draw();
  c1_up->cd();
  ///bkg///
  bkgstack_HN->Draw("hist");
  bkgstack_HN->SetMaximum(y_max[0]);
  bkgstack_HN->GetXaxis()->SetLabelSize(0);
  bkgstack_HN->GetYaxis()->SetTitle("Events / "+TString::Itoa(HN_dx,10)+" GeV");
  bkgstack_HN->GetYaxis()->SetLabelSize(0.05);
  bkgstack_HN->GetYaxis()->SetTitleSize(0.05);
  bkgstack_HN->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_HN[1]->SetMarkerStyle(3);
  hist_HN[1]->SetMarkerSize(1);
  hist_HN[1]->Draw("APsameE1");
  ///signal///
  hist_HN[0]->SetLineColor(kRed);
  hist_HN[0]->Draw("histsame");
  ///err///
  hist_HN_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_HN_err->SetFillStyle(3004);
  hist_HN_err->SetFillColor(kBlue);
  hist_HN_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c1_down->cd();
  TH1D *hist_HN_diff = (TH1D*)hist_HN[1]->Clone();
  hist_HN_diff->Divide(hist_HN_err);
  //hist_HN_diff->Add(hist_HN_err,-1);
  hist_HN_diff->SetMaximum(2);
  hist_HN_diff->SetMinimum(0);
  hist_HN_diff->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  hist_HN_diff->GetXaxis()->SetLabelSize(0.10);
  hist_HN_diff->GetXaxis()->SetTitleSize(0.10);
  hist_HN_diff->GetYaxis()->SetLabelSize(0.08);
  hist_HN_diff->SetYTitle("#frac{DATA}{MC}");
  hist_HN_diff->GetYaxis()->SetTitleSize(0.08);
  hist_HN_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_HN_diff->SetFillColorAlpha(45,0.35);
  hist_HN_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// W_real ///
  TCanvas *c2 = new TCanvas("c2", "", 800, 800);
  TPad *c2_up = new TPad("c2_up", "", 0, 0.25, 1, 1);
  c2_up->SetTopMargin( 0.05 ); c2_up->SetBottomMargin( 0.02 ); c2_up->SetRightMargin( 0.02 ); c2_up->SetLeftMargin( 0.1 );
  TPad *c2_down = new TPad("c2_down", "", 0, 0, 1, 0.25);
  c2_down->SetTopMargin( 0.03 ); c2_down->SetBottomMargin( 0.25 ); c2_down->SetRightMargin( 0.02 ); c2_down->SetLeftMargin( 0.1 ); c2_down->SetGridx(); c2_down->SetGridy();
  c2_up->Draw();
  c2_down->Draw();
  c2_up->cd();
  ///bkg///
  bkgstack_W_real->Draw("hist");
  bkgstack_W_real->SetMaximum(y_max[1]);
  bkgstack_W_real->GetXaxis()->SetLabelSize(0);
  bkgstack_W_real->GetYaxis()->SetTitle("Events / "+TString::Itoa(W_real_dx,10)+" GeV");
  bkgstack_W_real->GetYaxis()->SetLabelSize(0.05);
  bkgstack_W_real->GetYaxis()->SetTitleSize(0.05);
  bkgstack_W_real->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_W_real[1]->SetMarkerStyle(3);
  hist_W_real[1]->SetMarkerSize(1);
  hist_W_real[1]->Draw("APsameE1");
  ///signal///
  hist_W_real[0]->SetLineColor(kRed);
  hist_W_real[0]->Draw("histsame");
  ///err///
  hist_W_real_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_W_real_err->SetFillStyle(3004);
  hist_W_real_err->SetFillColor(kBlue);
  hist_W_real_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c2_down->cd();
  TH1D *hist_W_real_diff = (TH1D*)hist_W_real[1]->Clone();
  hist_W_real_diff->Divide(hist_W_real_err);
  //hist_W_real_diff->Add(hist_W_real_err,-1);
  hist_W_real_diff->SetMaximum(2);
  hist_W_real_diff->SetMinimum(0);
  hist_W_real_diff->GetXaxis()->SetTitle("m(#mu_{1}#mu_{2}#mu_{3}#nu) [GeV]");
  hist_W_real_diff->GetXaxis()->SetLabelSize(0.10);
  hist_W_real_diff->GetXaxis()->SetTitleSize(0.10);
  hist_W_real_diff->GetYaxis()->SetLabelSize(0.08);
  hist_W_real_diff->SetYTitle("#frac{DATA}{MC}");
  hist_W_real_diff->GetYaxis()->SetTitleSize(0.08);
  hist_W_real_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_W_real_diff->SetFillColorAlpha(45,0.35);
  hist_W_real_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// dR ///
  TCanvas *c3 = new TCanvas("c3", "", 800, 800);
  TPad *c3_up = new TPad("c3_up", "", 0, 0.25, 1, 1);
  c3_up->SetTopMargin( 0.05 ); c3_up->SetBottomMargin( 0.02 ); c3_up->SetRightMargin( 0.02 ); c3_up->SetLeftMargin( 0.1 );
  TPad *c3_down = new TPad("c3_down", "", 0, 0, 1, 0.25);
  c3_down->SetTopMargin( 0.03 ); c3_down->SetBottomMargin( 0.25 ); c3_down->SetRightMargin( 0.02 ); c3_down->SetLeftMargin( 0.1 ); c3_down->SetGridx(); c3_down->SetGridy();
  c3_up->Draw();
  c3_down->Draw();
  c3_up->cd();
  ///bkg///
  bkgstack_dR->Draw("hist");
  bkgstack_dR->SetMaximum(y_max[2]);
  bkgstack_dR->GetXaxis()->SetLabelSize(0);
  bkgstack_dR->GetYaxis()->SetTitle("Events / 0."+TString::Itoa(10*dR_dx,10));
  bkgstack_dR->GetYaxis()->SetLabelSize(0.05);
  bkgstack_dR->GetYaxis()->SetTitleSize(0.05);
  bkgstack_dR->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_dR[1]->SetMarkerStyle(3);
  hist_dR[1]->SetMarkerSize(1);
  hist_dR[1]->Draw("APsameE1");
  ///signal///
  hist_dR[0]->SetLineColor(kRed);
  hist_dR[0]->Draw("histsame");
  ///err///
  hist_dR_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_dR_err->SetFillStyle(3004);
  hist_dR_err->SetFillColor(kBlue);
  hist_dR_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c3_down->cd();
  TH1D *hist_dR_diff = (TH1D*)hist_dR[1]->Clone();
  hist_dR_diff->Divide(hist_dR_err);
  //hist_dR_diff->Add(hist_dR_err,-1);
  hist_dR_diff->SetMaximum(2);
  hist_dR_diff->SetMinimum(0);
  hist_dR_diff->GetXaxis()->SetTitle("#DeltaR(OS)_{Min}");
  hist_dR_diff->GetXaxis()->SetLabelSize(0.10);
  hist_dR_diff->GetXaxis()->SetTitleSize(0.10);
  hist_dR_diff->GetYaxis()->SetLabelSize(0.08);
  hist_dR_diff->SetYTitle("#frac{DATA}{MC}");
  hist_dR_diff->GetYaxis()->SetTitleSize(0.08);
  hist_dR_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_dR_diff->SetFillColorAlpha(45,0.35);
  hist_dR_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// gamma_star ///
  TCanvas *c4 = new TCanvas("c4", "", 800, 800);
  TPad *c4_up = new TPad("c4_up", "", 0, 0.25, 1, 1);
  c4_up->SetTopMargin( 0.05 ); c4_up->SetBottomMargin( 0.02 ); c4_up->SetRightMargin( 0.02 ); c4_up->SetLeftMargin( 0.1 );
  TPad *c4_down = new TPad("c4_down", "", 0, 0, 1, 0.25);
  c4_down->SetTopMargin( 0.03 ); c4_down->SetBottomMargin( 0.25 ); c4_down->SetRightMargin( 0.02 ); c4_down->SetLeftMargin( 0.1 ); c4_down->SetGridx(); c4_down->SetGridy();
  c4_up->Draw();
  c4_down->Draw();
  c4_up->cd();
  ///bkg///
  bkgstack_gamma_star->Draw("hist");
  bkgstack_gamma_star->SetMaximum(y_max[3]);
  bkgstack_gamma_star->GetXaxis()->SetLabelSize(0);
  bkgstack_gamma_star->GetYaxis()->SetTitle("Events / "+TString::Itoa(gamma_star_dx,10)+" GeV");
  bkgstack_gamma_star->GetYaxis()->SetLabelSize(0.05);
  bkgstack_gamma_star->GetYaxis()->SetTitleSize(0.05);
  bkgstack_gamma_star->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_gamma_star[1]->SetMarkerStyle(3);
  hist_gamma_star[1]->SetMarkerSize(1);
  hist_gamma_star[1]->Draw("APsameE1");
  ///signal///
  hist_gamma_star[0]->SetLineColor(kRed);
  hist_gamma_star[0]->Draw("histsame");
  ///err///
  hist_gamma_star_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_gamma_star_err->SetFillStyle(3004);
  hist_gamma_star_err->SetFillColor(kBlue);
  hist_gamma_star_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c4_down->cd();
  TH1D *hist_gamma_star_diff = (TH1D*)hist_gamma_star[1]->Clone();
  hist_gamma_star_diff->Divide(hist_gamma_star_err);
  //hist_gamma_star_diff->Add(hist_gamma_star_err,-1);
  hist_gamma_star_diff->SetMaximum(2);
  hist_gamma_star_diff->SetMinimum(0);
  hist_gamma_star_diff->GetXaxis()->SetTitle("m(#mu+#mu-) for Minimum #DeltaR [GeV]");
  hist_gamma_star_diff->GetXaxis()->SetLabelSize(0.10);
  hist_gamma_star_diff->GetXaxis()->SetTitleSize(0.10);
  hist_gamma_star_diff->GetYaxis()->SetLabelSize(0.08);
  hist_gamma_star_diff->SetYTitle("#frac{DATA}{MC}");
  hist_gamma_star_diff->GetYaxis()->SetTitleSize(0.08);
  hist_gamma_star_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_gamma_star_diff->SetFillColorAlpha(45,0.35);
  hist_gamma_star_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  /// z_candidate ///
  TCanvas *c5 = new TCanvas("c5", "", 800, 800);
  TPad *c5_up = new TPad("c5_up", "", 0, 0.25, 1, 1);
  c5_up->SetTopMargin( 0.05 ); c5_up->SetBottomMargin( 0.02 ); c5_up->SetRightMargin( 0.02 ); c5_up->SetLeftMargin( 0.1 );
  TPad *c5_down = new TPad("c5_down", "", 0, 0, 1, 0.25);
  c5_down->SetTopMargin( 0.03 ); c5_down->SetBottomMargin( 0.25 ); c5_down->SetRightMargin( 0.02 ); c5_down->SetLeftMargin( 0.1 ); c5_down->SetGridx(); c5_down->SetGridy();
  c5_up->Draw();
  c5_down->Draw();
  c5_up->cd();
  ///bkg///
  bkgstack_z_candidate->Draw("hist");
  bkgstack_z_candidate->SetMaximum(y_max[4]);
  bkgstack_z_candidate->GetXaxis()->SetLabelSize(0);
  bkgstack_z_candidate->GetYaxis()->SetTitle("Events / "+TString::Itoa(z_candidate_dx,10)+" GeV");
  bkgstack_z_candidate->GetYaxis()->SetLabelSize(0.05);
  bkgstack_z_candidate->GetYaxis()->SetTitleSize(0.05);
  bkgstack_z_candidate->GetYaxis()->SetTitleOffset(1.03);
  ///data///
  hist_z_candidate[1]->SetMarkerStyle(3);
  hist_z_candidate[1]->SetMarkerSize(1);
  hist_z_candidate[1]->Draw("APsameE1");
  ///signal///
  hist_z_candidate[0]->SetLineColor(kRed);
  hist_z_candidate[0]->Draw("histsame");
  ///err///
  cout << hist_z_candidate_err->GetXaxis()->GetNbins() << endl;
  hist_z_candidate_err->SetMarkerColorAlpha(kAzure-9, 0);
  hist_z_candidate_err->SetFillStyle(3004);
  hist_z_candidate_err->SetFillColor(kBlue);
  hist_z_candidate_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  ///MC-DATA///
  c5_down->cd();
  TH1D *hist_z_candidate_diff = (TH1D*)hist_z_candidate[1]->Clone();
  hist_z_candidate_diff->Divide(hist_z_candidate_err);
  //hist_z_candidate_diff->Add(hist_z_candidate_err,-1);
  hist_z_candidate_diff->SetMaximum(2);
  hist_z_candidate_diff->SetMinimum(0);
  hist_z_candidate_diff->GetXaxis()->SetTitle("m(#mu+#mu-) for Z-candidate [GeV]");
  hist_z_candidate_diff->GetXaxis()->SetLabelSize(0.10);
  hist_z_candidate_diff->GetXaxis()->SetTitleSize(0.10);
  hist_z_candidate_diff->GetYaxis()->SetLabelSize(0.08);
  hist_z_candidate_diff->SetYTitle("#frac{DATA}{MC}");
  hist_z_candidate_diff->GetYaxis()->SetTitleSize(0.08);
  hist_z_candidate_diff->GetYaxis()->SetTitleOffset(0.4);
  hist_z_candidate_diff->SetFillColorAlpha(45,0.35);
  hist_z_candidate_diff->Draw("PE1same");
  ///y=1 line///
  g1->Draw("same");
  
  
  if(isSave){
    c1->SaveAs(savepath+TString::Itoa(sample, 10)+"_HN.root");
    c2->SaveAs(savepath+TString::Itoa(sample, 10)+"_W_real.root");
    c3->SaveAs(savepath+TString::Itoa(sample, 10)+"_dR.root");
    c4->SaveAs(savepath+TString::Itoa(sample, 10)+"_gamma_star.root");
    c5->SaveAs(savepath+TString::Itoa(sample, 10)+"_z_candidate.root");
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
