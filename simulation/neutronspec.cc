void neutronspec() {
	ifstream in;
//	in.open("temp.txt");
	in.open("Ed5p36MeVGeant_all.txt");
	TH1F *hE = new TH1F("hE","",400,4,10);
	TH1F *hp = new TH1F("hp","",400,4,10);
	while(in.good()) {
		double En;
		TString instring;
		in>>instring>>En;
		if(instring=="NeutronE") hE->Fill(En);
		if(instring=="ProtonHit") hp->Fill(En);

	}
	hE->SetLineColor(kBlue);
	hE->Scale(hp->GetEntries()/hE->GetEntries());
	hE->Draw("HIST");
	hp->SetLineColor(kRed);
	hp->Draw("SAME");
}
