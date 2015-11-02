void cutflow(){
	TFile *rootfile[4];
	rootfile[0] = new TFile("ExampleAnalyzerElectronMuon_periodA_SKemu_dilep_5_3_14.root");
	rootfile[1] = new TFile("ExampleAnalyzerElectronMuon_periodB_SKemu_dilep_5_3_14.root");
	rootfile[2] = new TFile("ExampleAnalyzerElectronMuon_periodC_SKemu_dilep_5_3_14.root");
	rootfile[3] = new TFile("ExampleAnalyzerElectronMuon_periodD_SKemu_dilep_5_3_14.root");
	TH1F *cutflow[4];
	for(int i=0; i<4; i++){
		cutflow[i] = (TH1F*)rootfile[i]->Get("cutflow");
	}
	cutflow[0]->Add(cutflow[1]);
	cutflow[0]->Add(cutflow[2]);
	cutflow[0]->Add(cutflow[3]);
	cutflow[0]->Draw("hist");
	cutflow[0]->SetTitle("eemu_veto_loose");
}

