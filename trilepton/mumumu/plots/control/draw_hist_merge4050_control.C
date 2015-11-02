void draw_hist_merge4050_control(int muonidsel_int){
  
  gStyle->SetOptStat(0);
  
  TString var[5] = {"met", "lead_pt", "lead_eta", "n_jets", "n_bjets"};
  TString additional_cut[6] = {"before_dR", "dR", "dR_W", "dR_chi2", "dR_chi2_lambda", "dR_nbjets"};
  TString mass[3] = {"40", "50", "60"};
  TString muonidsel;
  TString additional_cut_latex[5] = { "#splitline{Before}{#splitline{Additional}{Cuts}}",
    "#splitline{#DeltaR > 0.5}{#splitline{ }{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{m(W_{RECO}) < 100 GeV}{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{#chi^{2} < 0.1}{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{#chi^{2} < 0.1}{#lambda = 0}}"};
  
  double y_max_loose[6*5] = {
    //  met     lead_pt   lead_eta  n_jets  n_bjets
    1000,   1000,     1000,     2000,   10000,      // before_dR
    800,    800,      800,      1200,   10000,      // dR
    100,    100,      100,      100,    1000,       // dR_W
    60,     80,       60,       60,     1000,       // dR_chi2
    60,     80,       60,       60,     1000,       // dR_chi2_lambda
    100,    100,      100,      150,    1000        // dR_nbjets
  };
  
  double y_max_tight[6*5] = {
    //  met     lead_pt   lead_eta  n_jets  n_bjets
    500,    600,      800,      600,    600,        // before_dR
    350,    600,      200,      200,    600,        // dR
    100,    100,      100,       60,     60,        // dR_W
    70,    100,       80,       40,     40,        // dR_chi2
    70,    100,       80,       40,     40,        // dR_chi2_lambda
    70,    100,       80,       40,     40         // dR_nbjets
  };
  
  
  double *y_max;
  if(muonidsel_int == 0){
    muonidsel = "loose";
    y_max = y_max_loose;
  }
  if(muonidsel_int == 1){
    muonidsel = "tight";
    y_max = y_max_tight;
  }
  
  
  for(int i=0;i<6;i++){ // cut
    double coupling_const;
    if(i==0 || i==1){
      coupling_const=0.001;
    }
    else{
      coupling_const=0.0001;
    }
    for(int j=0;j<5;j++){ // var
      TFile *file[3];
      TCanvas *can[3];
      file[0] = new TFile("./"+additional_cut[i]+"/"+muonidsel+"_40_"+var[j]+".root");
      file[1] = new TFile("./"+additional_cut[i]+"/"+muonidsel+"_50_"+var[j]+".root");
      file[2] = new TFile("./"+additional_cut[i]+"/"+muonidsel+"_60_"+var[j]+".root");
      can[0] = (TCanvas*)file[0]->Get("c"+TString::Itoa(j+1,10)); can[0]->Draw();
      can[1] = (TCanvas*)file[1]->Get("c"+TString::Itoa(j+1,10));
      can[2] = (TCanvas*)file[2]->Get("c"+TString::Itoa(j+1,10));
      TPad *pad[3];
      pad[0] = (TPad*)can[0]->GetPrimitive("c"+TString::Itoa(j+1,10)+"_up");
      pad[1] = (TPad*)can[1]->GetPrimitive("c"+TString::Itoa(j+1,10)+"_up");
      pad[2] = (TPad*)can[2]->GetPrimitive("c"+TString::Itoa(j+1,10)+"_up");
      pad[0]->cd();
      THStack *bkg_40 = (THStack*)pad[0]->GetPrimitive("bkgstack_"+var[j]);
      bkg_40->SetMaximum(y_max[5*i+j]);
      TH1D *hist_sig_50 = (TH1D*)pad[1]->GetPrimitive("hist_"+var[j]+"0");
      hist_sig_50->SetLineColor(kBlack);
      hist_sig_50->Draw("histsame");
      TH1D *hist_sig_60 = (TH1D*)pad[2]->GetPrimitive("hist_"+var[j]+"0");
      hist_sig_60->SetLineColor(kBlue);
      hist_sig_60->Draw("histsame");
      TLegend *lg = (TLegend*)pad[0]->GetPrimitive("TPave")->Clone();
      pad[0]->GetPrimitive("TPave")->Delete();
      lg->AddEntry(hist_sig_50, "HN50, |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
      lg->AddEntry(hist_sig_60, "HN60, |V_{N#mu}|^{2}=10^{"+TString::Itoa(TMath::Log10(coupling_const), 10)+"}", "l");
      lg->Draw();
      TLatex t;
      t.SetTextFont(43);
      t.SetTextColor(2);
      t.SetTextSize(20);
      double x_max = bkg_40->GetXaxis()->GetXmax(), x_min = bkg_40->GetXaxis()->GetXmin();
      if(j==4){
        bkg_40->SetMinimum(0.5);
        t.DrawLatex(x_min+0.72*(x_max-x_min), TMath::Power(10, 0.50*(TMath::Log10(y_max[5*i+j])-0.8)), additional_cut_latex[i]);
      }
      else{
        t.DrawLatex(x_min+0.72*(x_max-x_min), 0.50*y_max[5*i+j], additional_cut_latex[i]);
      }
      can[0]->SaveAs("./"+additional_cut[i]+"/"+muonidsel+"_"+TString::Itoa(405060,10)+"_"+var[j]+".pdf");
      can[0]->Close();
    }
  }
}