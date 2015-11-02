////// mumue channel //////

#include "/Users/jskim/Documents/CMS/header/Muon_charged_higgs.h"
#include "/Users/jskim/Documents/CMS/header/Electron.h"
#include "/Users/jskim/Documents/CMS/header/Event.h"
#include "/Users/jskim/Documents/CMS/header/Jet.h"
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

#define top_width 2.0
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

void charged_higgs_mumue(){
  TStopwatch totaltime;
  totaltime.Start();
  int sample, isPU, ischi2, isdR, isWreco, islambda, isSave=1;
  TString folder = "tri_2lightjet_1bjet";
  
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
    "TTG", "TTWW", "WWG", "WWW", "WWZ", "WZZ", "ZZZ", "ttZ", "HtoWW", "Wtollln_new"}; // others
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
  int N_entry, jet_index, solution_selection, smaller;
  struct Fit fitresult[2];
  Color_t fillcolors[N_chain] = {0, kRed-5, // signal and data. data random color
                                 kAzure-5, kAzure-5, kAzure-5,
                                 kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5, kGreen-5,
                                 kYellow-5,
                                 kRed-5,
                                 kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kGreen-10, kOrange-4};
 
  TH1D *hist_dark_photon[N_chain];
  TH1D *hist_dark_photon_err = new TH1D("hist_dark_photon_err", "", 120./1., 0, 120);
  THStack *bkgstack_dark_photon = new THStack("bkgstack_dark_photon", "");
  TH1D *hist_W_real[N_chain];
  TH1D *hist_W_real_err = new TH1D("hist_W_real_err", "", 500./10., 0, 500);
  THStack *bkgstack_W_real = new THStack("bkgstack_W_real", "");
  TH1D *hist_charged_higgs[N_chain];
  TH1D *hist_charged_higgs_err = new TH1D("hist_charged_higgs_err", "", 500./10., 0, 500);
  THStack *bkgstack_charged_higgs = new THStack("bkgstack_charged_higgs", "");
  TH1D *hist_top[N_chain];
  TH1D *hist_top_err = new TH1D("hist_top_err", "", 500./10., 0, 500);
  THStack *bkgstack_top = new THStack("bkgstack_top", "");
  
  TLorentzVector mu_4vec[2], el_4vec, dark_photon, W_real, charged_higgs, top, selection_nu[2];

  int FillOrder[N_chain] = { 0, 1, // signal, data
    14, 15, 16, 17, 18, 19, 20, 21, 22, // others
    12, // Wjets
    23, // Wtollln
    5, 6, 7, 8, 9, 10, 11, // VV
    13, // top
    2, 3, 4 }; // DY
  
  for(int l=1; l<N_chain; l++){
    int k = FillOrder[l];
    N_entry = 0;
    PUweight_sum=0;
 
    mu[k] = new Muon(chain_mu[k], 2);
    el[k] = new Electron(chain_el[k], 1);
    evt[k] = new Event(chain_event[k]);
    jet[k] = new Jet(chain_jet[k], 10);
    
    hist_dark_photon[k] = new TH1D("hist_dark_photon"+TString::Itoa(k,10), "", 120./1., 0, 120);
    hist_W_real[k] = new TH1D("hist_W_real"+TString::Itoa(k,10), "", 500./10., 0, 500);
    hist_charged_higgs[k] = new TH1D("hist_charged_higgs"+TString::Itoa(k,10), "", 500./10., 0, 500);
    hist_top[k] = new TH1D("hist_top"+TString::Itoa(k,10), "", 500./10., 0, 500);

    cout << endl << "Running Sample : " << chainlist[k] << endl;
    
    for(int i=0; i<chain_mu[k]->GetEntries(); i++){
      
      printeventnumber(i, chain_mu[k]->GetEntries());
      
      chain_mu[k]->GetEntry(i);
      chain_el[k]->GetEntry(i);
      chain_event[k]->GetEntry(i);
      chain_jet[k]->GetEntry(i);
      
      if( mu[k]->isHNLooseMuon(0) && mu[k]->isHNLooseMuon(1) && jet[k]->n_Jet<=10){
      //if( mu[k]->PassID("TightMuon", 0) && mu[k]->PassID("TightMuon", 1) && mu[k]->PassID("TightMuon", 2) ){
      //if( mu[k]->isHNTightMuon(0) && mu[k]->isHNTightMuon(1) && mu[k]->isHNTightMuon(2) ){
      //if( 1 ){ // NoHNLoose cut
        mu_4vec[0].SetPxPyPzE(mu[k]->Px[0],mu[k]->Py[0],mu[k]->Pz[0],mu[k]->E[0]);
        mu_4vec[1].SetPxPyPzE(mu[k]->Px[1],mu[k]->Py[1],mu[k]->Pz[1],mu[k]->E[1]);
        el_4vec.SetPxPyPzE(el[k]->Px[0],el[k]->Py[0],el[k]->Pz[0],el[k]->E[0]);
        
        MET = evt[k]->MET;
        METphi = evt[k]->METphi;
        
        /*
        if(ischi2){
          fit(&fitresult[0], 3, 1000, bjet+lep[0]+lep[1]+lep[2], MET, METphi, "m");
          fit(&fitresult[1], 3, 1000, bjet+lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 0 = minus, 1 = plus
          PutNuPz(&selection_nu[0], fitresult[0].Pz_fit);
          PutNuPz(&selection_nu[1], fitresult[1].Pz_fit);
        }
        else{
          PutNuPz(&selection_nu[0], solveqdeq(173.21, bjet+lep[0]+lep[1]+lep[2], MET, METphi, "m"));
          PutNuPz(&selection_nu[1], solveqdeq(173.21, bjet+lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus
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
        */
        
        int  n_light_jet = (jet[k]->n_Jet)-(evt[k]->n_bjet), n_bjet = (evt[k]->n_bjet);
        
        if(mu[k]->Charge[0] != mu[k]->Charge[1] && n_light_jet>=2 && n_bjet>=1){
        //if(mu[k]->Charge[0] != mu[k]->Charge[1] && n_light_jet>=2 && n_bjet>=1){
          
          dark_photon = mu_4vec[0]+mu_4vec[1];
          
          std::vector<TLorentzVector> bjet(n_bjet), light_jet(n_light_jet);

          /// find b-jet ///
          jet_index = 0;
          for(int j=0;j<jet[k]->n_Jet;j++){
            if(jet[k]->CombinedSecVertexBtag[j] > 0.679){
              //cout << j << '\t' << jet[k]->CombinedSecVertexBtag[j] << endl;
              bjet.at(jet_index).SetPxPyPzE(jet[k]->Px[j], jet[k]->Py[j], jet[k]->Pz[j], jet[k]->E[j]);
              jet_index++;
            }
          }
 
          /// find light jet ///
          jet_index = 0;
          for(int j=0;j<jet[k]->n_Jet;j++){
            if(jet[k]->CombinedSecVertexBtag[j] <= 0.679){
              light_jet.at(jet_index).SetPxPyPzE(jet[k]->Px[j], jet[k]->Py[j], jet[k]->Pz[j], jet[k]->E[j]);
              jet_index++;
            }
          }
          
          //W mass//
          double W_mass_temp=9999999;
          int good_W_jet_index[2];
          for(int j=0;j<n_light_jet-1;j++){
            for(int k=j+1;k<n_light_jet;k++){
              W_real = light_jet.at(j)+light_jet.at(k);
              if( fabs(W_real.M()-80.4) < fabs(W_mass_temp-80.4) ){
                W_mass_temp = W_real.M();
                good_W_jet_index[0]=j;
                good_W_jet_index[1]=k;
              }
            }
          }
          //W_real = light_jet.at(good_W_jet_index[0])+light_jet.at(good_W_jet_index[1]);
          W_real = light_jet.at(0)+light_jet.at(1);
          //W_real = light_jet.at(0)+bjet.at(0);
          
          charged_higgs = W_real + dark_photon;
          
          //top mass//
          double top_mass_temp=9999999;
          int good_bjet_index;
          for(int j=0;j<n_bjet;j++){
            top = charged_higgs+bjet.at(j);
            if( fabs(top.M()-173.21) < fabs(top_mass_temp-173.21) ){
              top_mass_temp = top.M();
              good_bjet_index=j;
            }
          }
          //top = charged_higgs+bjet.at(good_bjet_index);
          top = charged_higgs+bjet.at(0);
 
          if(!isPU) weight = evt[k]->weight;
          else      weight = evt[k]->weight * evt[k]->PU_reweight; // pileup weight
          
          /// apply weight event by event///
          //if(W_real.M() < 100.0){
          if(1){
            N_entry++;
            PUweight_sum+=evt[k]->PU_reweight;
            hist_dark_photon[k]->Fill(dark_photon.M(), weight);
            hist_W_real[k]->Fill(W_real.M(), weight);
            hist_charged_higgs[k]->Fill(charged_higgs.M(), weight);
            hist_top[k]->Fill(top.M(), weight);
            
          } // W < 100.0 GeV cut
          
          /*
          /// apply weight after filling histogram ///
          if(W_real.M() < 100.0){
          //if(1){
            N_entry++;
            PUweight_sum+=evt[k]->PU_reweight;
            hist_HN[k]->Fill(HN.M());
            hist_W_real[k]->Fill(W_real.M());
            hist_dR[k]->Fill( TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ));
            
            if(dark_photon.M() < 7.0)  hist_dark_photon[k]->Fill(dark_photon.M());
            else                      hist_dark_photon[k]->Fill(7.5);
          }
           */
        } // dR > 0.5 cut
      } // Loosemuon cut
    } // event loop

    cout << "weight = " << evt[k]->weight << endl
         << "entry = " << N_entry << endl
    << "PUweight = " << PUweight_sum/(double)N_entry << endl;

    /// apply weight after filling histogram ///
    //hist_HN[k]->Scale(weight);
    //hist_W_real[k]->Scale(weight);
    //hist_dR[k]->Scale(weight);
    //hist_dark_photon[k]->Scale(weight);

    /*
    if(k==4){ // Add two DY process
      hist_HN[4]->Add(hist_HN[3]);
      hist_W_real[4]->Add(hist_W_real[3]);
    }
    */
    
    hist_dark_photon[k]->SetFillColor(fillcolors[k]);
    hist_W_real[k]->SetFillColor(fillcolors[k]);
    hist_charged_higgs[k]->SetFillColor(fillcolors[k]);
    hist_top[k]->SetFillColor(fillcolors[k]);
 
    /// merging same-type histograms
    if(k>=15 && k<=22){
      hist_dark_photon[14]->Add(hist_dark_photon[k]);
      hist_W_real[14]->Add(hist_W_real[k]);
      hist_charged_higgs[14]->Add(hist_charged_higgs[k]);
      hist_top[14]->Add(hist_top[k]);
    }
    if(k>=6 && k<=11){
      hist_dark_photon[5]->Add(hist_dark_photon[k]);
      hist_W_real[5]->Add(hist_W_real[k]);
      hist_charged_higgs[5]->Add(hist_charged_higgs[k]);
      hist_top[5]->Add(hist_top[k]);
    }
    if(k==3 || k==4){
      hist_dark_photon[2]->Add(hist_dark_photon[k]);
      hist_W_real[2]->Add(hist_W_real[k]);
      hist_charged_higgs[2]->Add(hist_charged_higgs[k]);
      hist_top[2]->Add(hist_top[k]);
    }
    if(k==14 || k==23 || k==5 || k==12 || k==13 || k==2){
      bkgstack_dark_photon->Add(hist_dark_photon[k]);
      bkgstack_W_real->Add(hist_W_real[k]);
      bkgstack_charged_higgs->Add(hist_charged_higgs[k]);
      bkgstack_top->Add(hist_top[k]);
    }
    if(k>=2){
      hist_dark_photon_err->Add(hist_dark_photon[k]);
      hist_W_real_err->Add(hist_W_real[k]);
      hist_charged_higgs_err->Add(hist_charged_higgs[k]);
      hist_top_err->Add(hist_top[k]);
    }

  
  } // chain loop
  /*
  /// signal scale factor ///
  double k_factor = 1.34, coupling_const;
  if(sample==40)  coupling_const=0.0001;
  if(sample==50)  coupling_const=0.0001;
  if(sample==60)  coupling_const=0.0001;
  if(sample==100)  coupling_const=0.01;
  hist_dark_photon[0]->Scale(k_factor*coupling_const);
  hist_W_real[0]->Scale(k_factor*coupling_const);
  hist_charged_higgs[0]->Scale(k_factor*coupling_const);
  hist_top[0]->Scale(k_factor*coupling_const);
  */
  
  Double_t totalRuntime = totaltime.CpuTime();
  cout << "\n----------Total Runtime: " << totalRuntime << " seconds----------\n" << endl;
  
  gStyle->SetOptStat(0);
  
  /// top ///
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->Draw();
  ///bkg///
  bkgstack_top->Draw("hist");
  bkgstack_top->SetMaximum(40); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_top->SetTitle("top");
  //bkgstack_top->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  bkgstack_top->GetYaxis()->SetTitle("Event / 10 GeV");
  ///data///
  hist_top[1]->SetMarkerStyle(3);
  hist_top[1]->SetMarkerSize(1);
  hist_top[1]->Draw("APsameE1");
  ///signal///
  //hist_top[0]->SetLineColor(kRed);
  //hist_top[0]->Draw("histsame");
  ///err///
  hist_top_err->SetFillStyle(3004);
  hist_top_err->SetFillColor(kBlue);
  hist_top_err->Draw("sameE2");
  ///legend///
  TLegend *lg1 = new TLegend(0.65, 0.6, 0.9, 0.9);
  lg1->AddEntry(hist_top[1], "data", "p");
  //lg1->AddEntry(hist_top[0], "top"+TString::Itoa(sample, 10)+", |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
  lg1->AddEntry(hist_top[2], "DY", "f");
  lg1->AddEntry(hist_top[13], "top", "f");
  lg1->AddEntry(hist_top[5], "VV", "f");
  lg1->AddEntry(hist_top[23], "Wtollln", "f");
  lg1->AddEntry(hist_top[12], "Wjets", "f");
  lg1->AddEntry(hist_top[14], "other", "f");
  lg1->Draw();
  if(isSave){
    c1->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/top_mumue.root");
    c1->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/top_mumue.pdf");
  }
  
  /// charged_higgs ///
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  c2->Draw();
  ///bkg///
  bkgstack_charged_higgs->Draw("hist");
  bkgstack_charged_higgs->SetMaximum(50); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_charged_higgs->SetTitle("charged_higgs");
  //bkgstack_charged_higgs->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  bkgstack_charged_higgs->GetYaxis()->SetTitle("Event / 10 GeV");
  ///data///
  hist_charged_higgs[1]->SetMarkerStyle(3);
  hist_charged_higgs[1]->SetMarkerSize(1);
  hist_charged_higgs[1]->Draw("APsameE1");
  ///signal///
  //hist_charged_higgs[0]->SetLineColor(kRed);
  //hist_charged_higgs[0]->Draw("histsame");
  ///err///
  hist_charged_higgs_err->SetFillStyle(3004);
  hist_charged_higgs_err->SetFillColor(kBlue);
  hist_charged_higgs_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  if(isSave){
    c2->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/charged_higgs_mumue.root");
    c2->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/charged_higgs_mumue.pdf");
  }
  
  /// dark_photon ///
  TCanvas *c3 = new TCanvas("c3", "", 800, 600);
  c3->Draw();
  ///bkg///
  bkgstack_dark_photon->Draw("hist");
  bkgstack_dark_photon->SetMaximum(30); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_dark_photon->SetTitle("dark_photon");
  //bkgstack_dark_photon->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  bkgstack_dark_photon->GetYaxis()->SetTitle("Event / 1 GeV");
  ///data///
  hist_dark_photon[1]->SetMarkerStyle(3);
  hist_dark_photon[1]->SetMarkerSize(1);
  hist_dark_photon[1]->Draw("APsameE1");
  ///signal///
  //hist_dark_photon[0]->SetLineColor(kRed);
  //hist_dark_photon[0]->Draw("histsame");
  ///err///
  hist_dark_photon_err->SetFillStyle(3004);
  hist_dark_photon_err->SetFillColor(kBlue);
  hist_dark_photon_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  if(isSave){
    c3->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/dark_photon_mumue.root");
    c3->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/dark_photon_mumue.pdf");
  }
  
  /// W_real ///
  TCanvas *c4 = new TCanvas("c4", "", 800, 600);
  c4->Draw();
  ///bkg///
  bkgstack_W_real->Draw("hist");
  bkgstack_W_real->SetMaximum(60); // 100 for dR_W, 70 for chi2_lambda
  bkgstack_W_real->SetTitle("W_real");
  //bkgstack_W_real->GetXaxis()->SetTitle("m(#mu_{2}#mu_{3}#nu) [GeV]");
  bkgstack_W_real->GetYaxis()->SetTitle("Event / 10 GeV");
  ///data///
  hist_W_real[1]->SetMarkerStyle(3);
  hist_W_real[1]->SetMarkerSize(1);
  hist_W_real[1]->Draw("APsameE1");
  ///signal///
  //hist_W_real[0]->SetLineColor(kRed);
  //hist_W_real[0]->Draw("histsame");
  ///err///
  hist_W_real_err->SetFillStyle(3004);
  hist_W_real_err->SetFillColor(kBlue);
  hist_W_real_err->Draw("sameE2");
  ///legend///
  lg1->Draw();
  if(isSave){
    c4->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/W_real_mumue.root");
    c4->SaveAs("/Users/jskim/Documents/CMS/charged_higgs/"+folder+"/W_real_mumue.pdf");
  }
}



double chisquare(double M_fit, double lambda){
  return TMath::Power((M_fit-80.4)/top_width, 2) + lambda*lambda;
}

void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, TString pm){
  double dlambda = range_lambda/(n_lambda-1), lambda, x, x_temp, MET_temp, Pz_temp; // step size
  TLorentzVector nu, W_real;
  
  lambda = 0;
  MET_temp = MET_recon;
  nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
  Pz_temp = solveqdeq(&(fitresult->isSamePz), 173.21, l1l2l3, MET_temp, METphi_recon, pm);
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
      Pz_temp = solveqdeq(&(fitresult->isSamePz), 173.21, l1l2l3, MET_temp, METphi_recon, pm);
      PutNuPz(&nu, Pz_temp);
      W_real = l1l2l3 + nu;
      x_temp = chisquare(W_real.M(), lambda);
      if(W_real.M()==173.21) break; // lambda = 0 solvable
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
