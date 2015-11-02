TString frac_2p2(double x);

void draw_sig_recon(int muonidsel_int){
  
  gStyle->SetOptStat(0);
  
  TString var = "sig_recon";
  TString additional_cut[5] = {"before_dR", "dR", "dR_W", "dR_chi2", "dR_chi2_lambda"};
  TString mass[3] = {"40", "50", "60"};
  TString muonidsel;
  TString additional_cut_latex[5] = { "#splitline{Before}{#splitline{Additional}{Cuts}}",
    "#splitline{#DeltaR > 0.5}{#splitline{ }{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{m(W_{RECO}) < 100 GeV}{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{#chi^{2} < 0.1}{ }}",
    "#splitline{#DeltaR > 0.5}{#splitline{#chi^{2} < 0.1}{#lambda = 0}}"};
  Color_t linecolor[3] = {kBlack, kRed, kBlue};
  
  if(muonidsel_int == 0){
    muonidsel = "loose";
  }
  if(muonidsel_int == 1){
    muonidsel = "tight";
  }
    
    
  for(int i=0;i<5;i++){ // cut
    TFile *file[3];
    TCanvas *can[3];
    TH1D *hist[3];
    TLegend lg(0.9, 0.7, 1, 0.9);
    TCanvas c("c","",800,600);
    c.Draw();
    for(int j=0;j<3;j++){ // mass
      file[j] = new TFile("./"+additional_cut[i]+"/"+muonidsel+"_"+mass[j]+"_"+var+".root");
      can[j] = (TCanvas*)file[j]->Get("c6");
      c.cd();
      hist[j] = (TH1D*)can[j]->GetPrimitive("hist_"+var);
      hist[j]->SetLineColor(linecolor[j]);
      hist[j]->Scale(1./hist[j]->Integral());
      hist[j]->Draw("histsame");
      lg.AddEntry(hist[j], "HN"+mass[j], "l");
      if(j==2)  lg.Draw();
      if(j==0){
        hist[0]->SetMaximum(0.1);
        hist[0]->SetXTitle("m(HN_{RECO}) [GeV]");
      }
    }
    TLatex t, statinfo;
    t.SetTextFont(43);
   	t.SetTextColor(2);
    t.SetTextSize(20);
    t.DrawLatex(100, 0.05, additional_cut_latex[i]);
    statinfo.SetTextFont(43);
    statinfo.SetTextSize(25);
    statinfo.DrawLatex(80, 0.085, "#splitline{           Mean  Width}{#splitline{HN40  "+frac_2p2(hist[0]->GetMean())+"  "+frac_2p2(hist[0]->GetRMS())+"}{#splitline{HN50  "+frac_2p2(hist[1]->GetMean())+"  "+frac_2p2(hist[1]->GetRMS())+"}{HN60  "+frac_2p2(hist[2]->GetMean())+"  "+frac_2p2(hist[2]->GetRMS())+"}}}");
    c.SaveAs("./"+additional_cut[i]+"/"+muonidsel+"_"+var+".pdf");
  }
}

TString frac_2p2(double x){
  int A = (int)x, B = 100*x-100*A;
  if(x-(A+0.01*B) >= 0.05)  B++;
  if(B>=10) return TString::Itoa(A,10)+"."+TString::Itoa(B,10);
  else  return TString::Itoa(A,10)+".0"+TString::Itoa(B,10);
}
