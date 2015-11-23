#include "/Users/jskim/Documents/CMS/header/Muon.h"
#include "/Users/jskim/Documents/CMS/header/Event.h"
#include "/Users/jskim/Documents/CMS/header/Gen.h"
#include "/Users/jskim/Documents/CMS/header/trilepton.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <TRandom.h>

#define W_width 2.1
#define MET_width 7

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
void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, double W_mass, TString pm);


/////////////////////
/// main function ///
/////////////////////

void signal_recon(){
  int ischi2 = 0, additional_cut, sample_mass[3] = {40, 50, 60}, larger, smaller, l3nu_deltaR_larger, l3nu_deltaR_smaller, better, solution_selection, n_pass_cut[3]={0,0,0};
  double MET, METphi, chi2_sol_sel;
  struct Fit fitresult[2];
  Color_t linecolor[3] = {kBlack, kBlue, kRed};
  bool selection;
  
  TChain *chain_mu[3];
  TChain *chain_event[3];
  TChain *chain_gen[3];
  Muon *mu[3];
  Event *evt[3];
  Gen *gen[3];
  TLorentzVector lep[3], nu, HN, W_real, gen_lep[3], gen_nu, gen_HN, gen_W;
  TLorentzVector selection_nu[2], selection_HN[2];
  TH1D *hist_HN[3];
  TH2D *hist_Pz_fit_vs_true[3];
  TH1D *hist_larger_smaller[3];
  TH1D *hist_plus_minus[3];
  TH1D *hist_l3nu_deltaR_larger_smaller[3];
  TH1D *hist_sol_diff[3];
  TH1D *hist_new_method[3];
	TH1D *hist_gen_HN_mass[3];
  TH1D *hist_chi2[3];
  TH1D *hist_lambda[3];
  
  TH1D *hist_MET[3];
  TH1D *hist_gen_nu_pt[3];
  TH1D *hist_Pz_sol[3];
  TH1D *hist_gen_nu_pz[3];
  
  //gStyle->SetOptStat(0);
  
  cout
  << "===============================================" << endl
  << "Additional Cut?" << endl
  << "  0 : No Additional Cut" << endl
  << "  1 : delta R > 0.5" << endl
  << "  2 : delta R > 0.5 & W_reco < 100 GeV" << endl
  << "  3 : delta R > 0.5 & chi2 < 0.1" << endl
  << "  4 : delta R > 0.5 & chi2 < 0.1 & lambda = 0" << endl
  << "===============================================" << endl
  << "Cut : ";
  cin >> additional_cut;
  cout << endl;
  
  if(additional_cut == 3 || additional_cut == 4) ischi2 = 1;
  
  for(int sample_it=0; sample_it<3; sample_it++){
    chi2_sol_sel = 0;
    
    chain_mu[sample_it] = new TChain("Muon");
    chain_mu[sample_it]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample_mass[sample_it], 10)+"_new_5_3_14.root");
    chain_event[sample_it] = new TChain("Event");
    chain_event[sample_it]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample_mass[sample_it], 10)+"_new_5_3_14.root");
    chain_gen[sample_it] = new TChain("Gen");
    chain_gen[sample_it]->Add("./files/signal/trilepton_mumumu_SKHN"+TString::Itoa(sample_mass[sample_it], 10)+"_new_5_3_14.root");
    
    mu[sample_it] = new Muon(chain_mu[sample_it], 3);
    evt[sample_it] = new Event(chain_event[sample_it]);
    gen[sample_it] = new Gen(chain_gen[sample_it]);
    hist_HN[sample_it] = new TH1D("", "HN_{reco} invariant mass : #chi^{2} < 0.1 / #DeltaR > 0.5;m(#mu#mu#nu) [GeV];Normalized Entries", 100./1., 0., 100.);
    hist_Pz_fit_vs_true[sample_it] = new TH2D("", "#nu^{fit}_{L} vs #nu^{true}_{L} w/ HN"+TString::Itoa(sample_mass[sample_it], 10), 25./1., -50., 50., 25./1., -50., 50.);
    hist_larger_smaller[sample_it] = new TH1D("", "smaller vs larger", 2, 0, 2);
    hist_plus_minus[sample_it] = new TH1D("", "minus vs plus", 2, 0, 2);
    hist_l3nu_deltaR_larger_smaller[sample_it] = new TH1D("", "l3nu delta R smaller vs larger", 2, 0, 2);
    hist_sol_diff[sample_it] = new TH1D("", "| P_{z}+ #minus P_{z}- |", 100./1., 0., 100.);
    hist_gen_HN_mass[sample_it] = new TH1D("", "HN_{gen} invariant mass;m(#mu#mu#nu) [GeV];Normalized Entries", 30./0.1, 35., 65.);
    hist_chi2[sample_it] = new TH1D("", "#chi^{2}_{fit} w/ HN"+TString::Itoa(sample_mass[sample_it], 10), 20, 0., 20.);
    hist_lambda[sample_it] = new TH1D("", "#lambda_{fit} w/ HN"+TString::Itoa(sample_mass[sample_it], 10), 6./0.1, -3, 3);
    
    hist_MET[sample_it] = new TH1D("", "hist_MET", 100, 0, 100);
    hist_gen_nu_pt[sample_it] = new TH1D("", "hist_gen_nu_pt", 100, 0, 100);
    hist_Pz_sol[sample_it] = new TH1D("", "hist_Pz_sol", 100, 0, 100);
    hist_gen_nu_pz[sample_it] = new TH1D("", "hist_gen_nu_pz", 100, 0, 100);
    
    int lep_order[3] = {0, 1, 2};
    if(sample_mass[sample_it] >= 60){
      lep_order[0] = 2;
      lep_order[2] = 0;
    }
    
    for(int i=0; i<chain_mu[sample_it]->GetEntries(); i++){
      chain_mu[sample_it]->GetEntry(i);
      chain_event[sample_it]->GetEntry(i);
      chain_gen[sample_it]->GetEntry(i);
      
      for(int j=0; j<3; j++){
        lep[lep_order[j]].SetPxPyPzE(mu[sample_it]->Px[j],mu[sample_it]->Py[j],mu[sample_it]->Pz[j],mu[sample_it]->E[j]); // for reco lep, lep[0].pt() > lep[2].pt()
        gen_lep[j].SetPxPyPzE(gen[sample_it]->mu_Px[j],gen[sample_it]->mu_Py[j],gen[sample_it]->mu_Pz[j],gen[sample_it]->mu_E[j]); // for gen lep, l1 l2 l3 are in time-ordered in production
      }
      
      //Cut : Three HNLosse//
      selection = mu[sample_it]->isHNLooseMuon(0, 15) && mu[sample_it]->isHNLooseMuon(1, 10) && mu[sample_it]->isHNLooseMuon(2, 10);
      //Cut : dR > 0.5 //
      if(additional_cut == 1 || additional_cut == 2 || additional_cut == 3 || additional_cut == 4){
        selection *= TMath::Min( lep[0].DeltaR(lep[1]) , lep[2].DeltaR(lep[1]) ) > 0.5;
      }
      
      if(selection){
        gen_nu.SetPxPyPzE(gen[sample_it]->nu_Px, gen[sample_it]->nu_Py, gen[sample_it]->nu_Pz, gen[sample_it]->nu_E);
        gen_HN.SetPxPyPzE(gen[sample_it]->HN_Px, gen[sample_it]->HN_Py, gen[sample_it]->HN_Pz, gen[sample_it]->HN_E);
        gen_W.SetPxPyPzE(gen[sample_it]->W_Px, gen[sample_it]->W_Py, gen[sample_it]->W_Pz, gen[sample_it]->W_E);
        hist_gen_HN_mass[sample_it]->Fill(gen_HN.M());
        MET = evt[sample_it]->MET;
        METphi = evt[sample_it]->METphi;
        nu.SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        for(int j=0; j<2; j++){
          selection_nu[j].SetPxPyPzE(MET*TMath::Cos(METphi), MET*TMath::Sin(METphi), 0, 0);
        }
        
        if(ischi2){
          fit(&fitresult[0], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "m");
          fit(&fitresult[1], 3, 1000, lep[0]+lep[1]+lep[2], MET, METphi, "p"); // 0 = minus, 1 = plus
          //fit(&fitresult[0], 3, 1000, lep[0]+lep[1]+lep[2], gen_nu.Pt(), gen_nu.Phi(), 80.4, "m");
          //fit(&fitresult[1], 3, 1000, lep[0]+lep[1]+lep[2], gen_nu.Pt(), gen_nu.Phi(), 80.4, "p"); // 0 = minus, 1 = plus
          PutNuPz(&selection_nu[0], fitresult[0].Pz_fit);
          PutNuPz(&selection_nu[1], fitresult[1].Pz_fit);
        }
        else{
          PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "m"));
          PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], MET, METphi, "p")); // 0 = minus, 1 = plus
          //PutNuPz(&selection_nu[0], solveqdeq(80.4, lep[0]+lep[1]+lep[2], gen_nu.Pt(), gen_nu.Phi(), "m"));
          //PutNuPz(&selection_nu[1], solveqdeq(80.4, lep[0]+lep[1]+lep[2], gen_nu.Pt(), gen_nu.Phi(), "p")); // 0 = minus, 1 = plus
        }

        if(selection_nu[0].Pz() == selection_nu[1].Pz()){
          solution_selection = 0; // 0, 1 상관 없으므로.
        }
        else{         // different Pz에 대해서만 비교.
          // find the one with larger maginitude
          if(fabs(selection_nu[0].Pz()) > fabs(selection_nu[1].Pz())){
            larger = 0;
            smaller = 1;
          }
          else{
            larger = 1;
            smaller = 0;
          }
          // find the one with larger l3nu angle //
          if(selection_nu[0].DeltaR(lep[2]) > selection_nu[1].DeltaR(lep[2])){
            l3nu_deltaR_larger = 0;
            l3nu_deltaR_smaller = 1;
          }
          else{
            l3nu_deltaR_larger = 1;
            l3nu_deltaR_smaller = 0;
          }
          // check the better solution
          // in the histogram, 0 : smaller, 1 : larger
          if(fabs(gen_nu.Pz()-selection_nu[smaller].Pz()) < fabs(gen_nu.Pz()-selection_nu[larger].Pz())){
            better = smaller;
            hist_larger_smaller[sample_it]->Fill(0);
          }
          else{
            better = larger;
            hist_larger_smaller[sample_it]->Fill(1);
          }
          // in the histogram, 0 : minus, 1 : plus
          if(better == 0) hist_plus_minus[sample_it]->Fill(0);
          else            hist_plus_minus[sample_it]->Fill(1);
          // in the histogram, 0 : smaller, 1 : larger
          if(better == l3nu_deltaR_larger) hist_l3nu_deltaR_larger_smaller[sample_it]->Fill(1);
          else                             hist_l3nu_deltaR_larger_smaller[sample_it]->Fill(0);

          hist_sol_diff[sample_it]->Fill(fabs(selection_nu[0].Pz()-selection_nu[1].Pz()));
          
          solution_selection = smaller;
        }
        
        PutNuPz(&nu, selection_nu[solution_selection].Pz());
        HN = lep[1] + lep[2] + nu;
        W_real = HN + lep[0];
        
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
        
        if(selection){
          
          hist_HN[sample_it]->Fill(HN.M());
          hist_Pz_fit_vs_true[sample_it]->Fill(gen_nu.Pz(), selection_nu[solution_selection].Pz());
          if(!(selection_nu[0].Pz() == selection_nu[1].Pz())){
            n_pass_cut[sample_it]++;
            hist_chi2[sample_it]->Fill(TMath::Power( ( gen_nu.Pz() - selection_nu[solution_selection].Pz() )/gen_nu.Pz(), 2 ));
          }
          //hist_chi2[sample_it]->Fill(fitresult[solution_selection].chi2);
          hist_lambda[sample_it]->Fill(fitresult[solution_selection].lambda_fit);
          hist_MET[sample_it]->Fill(MET);
          hist_gen_nu_pt[sample_it]->Fill(gen_nu.Pt());
          hist_Pz_sol[sample_it]->Fill(nu.Pz());
          hist_gen_nu_pz[sample_it]->Fill(gen_nu.Pz());
          
        } // chi2 cut loop
      } // triloose loop
    } // event loop
    hist_larger_smaller[sample_it]->Scale(100./hist_larger_smaller[sample_it]->Integral());
    hist_plus_minus[sample_it]->Scale(100./hist_plus_minus[sample_it]->Integral());
    hist_l3nu_deltaR_larger_smaller[sample_it]->Scale(100./hist_l3nu_deltaR_larger_smaller[sample_it]->Integral());
    
    cout << "HN" << sample_mass[sample_it] << " : " << endl
    << '\t' << "pass cut : " << n_pass_cut[sample_it] << endl
    << '\t' << "mean : " << hist_HN[sample_it]->GetMean(1) << endl
    << '\t' << "rms  : " << hist_HN[sample_it]->GetRMS(1) << endl
    << '\t' << "plus : " << hist_plus_minus[sample_it]->GetBinContent(2) << endl
    << '\t' << "smaller : " << hist_larger_smaller[sample_it]->GetBinContent(1) << endl
    << '\t' << "l3nu deltaR larger : " << hist_l3nu_deltaR_larger_smaller[sample_it]->GetBinContent(2) << endl
    << '\t' << "chi2 for sol sel : " << hist_chi2[sample_it]->GetMean(1) << endl;
  } // sample_it loop
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->cd();
  for(int sample_it=0; sample_it<3; sample_it++){
    hist_HN[sample_it]->SetLineColor(linecolor[sample_it]);
    hist_HN[sample_it]->Scale(100./hist_HN[sample_it]->Integral());
    hist_HN[sample_it]->Draw("histsame");
  }
  hist_HN[0]->SetMaximum(15);
  TLegend *lg1 = new TLegend(0.9, 0.7, 1, 0.9);
  for(int i=0; i<3; i++)	lg1->AddEntry(hist_HN[i], "HN"+TString::Itoa(sample_mass[i], 10), "l");
  lg1->Draw();
  
  
  
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  c2->cd();
  hist_Pz_fit_vs_true[0]->Draw("scat");
  //hist_Pz_fit_vs_true[0]->SetTitle("#nu^{fit}_{L} vs #nu^{true}_{L} w/ HN40 w/ Pz_{better}");
  hist_Pz_fit_vs_true[0]->SetXTitle("Pz_{true}");
  hist_Pz_fit_vs_true[0]->SetYTitle("Pz_{fit}");
  

  
  TCanvas *c41 = new TCanvas("c41", "", 800, 600);
  c41->cd();
  hist_chi2[0]->Draw("hist");
  //c41->SaveAs("/Users/jskim/Desktop/plots/chi2_noaddcut_HN40.pdf");
  //TCanvas *c42 = new TCanvas("c42", "", 800, 600);
  //c42->cd();
  //hist_lambda[0]->Draw("hist");
  //c42->SaveAs("/Users/jskim/Desktop/plots/lambda_noaddcut_HN40.pdf");
  
  TCanvas *c51 = new TCanvas("c51", "", 800, 600);
  c51->cd();
  hist_chi2[1]->Draw("hist");
  //c51->SaveAs("/Users/jskim/Desktop/plots/chi2_noaddcut_HN50.pdf");
  //TCanvas *c52 = new TCanvas("c52", "", 800, 600);
  //c52->cd();
  //hist_lambda[1]->Draw("hist");
  //c52->SaveAs("/Users/jskim/Desktop/plots/lambda_noaddcut_HN50.pdf");
  
  TCanvas *c61 = new TCanvas("c61", "", 800, 600);
  c61->cd();
  hist_chi2[2]->Draw("hist");
  //c61->SaveAs("/Users/jskim/Desktop/plots/chi2_noaddcut_HN60.pdf");
  //TCanvas *c62 = new TCanvas("c62", "", 800, 600);
  //c62->cd();
  //hist_lambda[2]->Draw("hist");
  //c62->SaveAs("/Users/jskim/Desktop/plots/lambda_noaddcut_HN60.pdf");
  
  /*
  cout << "======GEN======" << endl;
  TCanvas *c5 = new TCanvas("c5", "", 800, 600);
  c5->cd();
  for(int sample_it=0; sample_it<3; sample_it++){
    hist_gen_HN_mass[sample_it]->SetLineColor(linecolor[sample_it]);
    hist_gen_HN_mass[sample_it]->Scale(100./hist_gen_HN_mass[sample_it]->Integral());
    hist_gen_HN_mass[sample_it]->Draw("histsame");
    cout << '\t' << "mean : " << hist_gen_HN_mass[sample_it]->GetMean(1) << endl;
    cout << '\t' << "rms  : " << hist_gen_HN_mass[sample_it]->GetRMS(1) << endl;
  }
  hist_gen_HN_mass[0]->SetMaximum(100);
  lg1->Draw();
  */

  
  
  TCanvas *c6 = new TCanvas("c6", "", 800, 600);
  hist_MET[0]->DrawNormalized("hist", 1);
  hist_gen_nu_pt[0]->DrawNormalized("psame", 1);
  c6->SaveAs("Pt.png");
  
  TCanvas *c7 = new TCanvas("c7", "", 800, 600);
  hist_Pz_sol[0]->DrawNormalized("hist", 1);
  hist_gen_nu_pz[0]->DrawNormalized("psame", 1);
  c7->SaveAs("Pz.png");
  
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

void fit(struct Fit *fitresult, double range_lambda, double n_lambda, TLorentzVector const l1l2l3, double MET_recon, double METphi_recon, double W_mass, TString pm){
  double dlambda = range_lambda/(n_lambda-1), lambda, x, x_temp, MET_temp, Pz_temp; // step size
  TLorentzVector nu, W_real;
  
  lambda = 0;
  MET_temp = MET_recon;
  nu.SetPxPyPzE(MET_temp*TMath::Cos(METphi_recon), MET_temp*TMath::Sin(METphi_recon), 0, 0);
  Pz_temp = solveqdeq(&(fitresult->isSamePz), W_mass, l1l2l3, MET_temp, METphi_recon, pm);
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
      Pz_temp = solveqdeq(&(fitresult->isSamePz), W_mass, l1l2l3, MET_temp, METphi_recon, pm);
      PutNuPz(&nu, Pz_temp);
      W_real = l1l2l3 + nu;
      x_temp = chisquare(W_real.M(), lambda);
      
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
