void draw_hist_control(int muonidsel_int){
  
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
    for(int j=0;j<5;j++){ // var
      for(int k=0;k<3;k++){ // mass
        TFile file("./"+additional_cut[i]+"/"+muonidsel+"_"+mass[k]+"_"+var[j]+".root");
        TCanvas *can = (TCanvas*)file.Get("c"+TString::Itoa(j+1,10));
        TPad *pad = (TPad*)can->GetPrimitive("c"+TString::Itoa(j+1,10)+"_up");
        THStack *bkg = (THStack*)pad->GetPrimitive("bkgstack_"+var[j]);
        bkg->SetMaximum(y_max[5*i+j]);
        TLatex t;
        t.SetTextFont(43);
        t.SetTextColor(2);
        t.SetTextSize(20);
        double x_max = bkg->GetXaxis()->GetXmax(), x_min = bkg->GetXaxis()->GetXmin();
        if(j==4){
          bkg->SetMinimum(0.5);
          t.DrawLatex(x_min+0.72*(x_max-x_min), TMath::Power(10, 0.50*(TMath::Log10(y_max[5*i+j])-0.8)), additional_cut_latex[i]);
        }
        else{
          t.DrawLatex(x_min+0.72*(x_max-x_min), 0.50*y_max[5*i+j], additional_cut_latex[i]);
        }
        
        can->SaveAs("./"+additional_cut[i]+"/"+muonidsel+"_"+mass[k]+"_"+var[j]+".pdf");
        can->Close();
      }
    }
  }
}