void plot_track() {
	TFile *f = TFile::Open("output.root");
	TTree *tree = (TTree*)f->Get("simData");
	std::vector<double> *pE = 0;
	std::vector<double> *px = 0;
	std::vector<double> *py = 0;
	std::vector<double> *pz = 0;
	tree->SetBranchAddress("pE",&pE);
	tree->SetBranchAddress("px",&px);
	tree->SetBranchAddress("py",&py);
	tree->SetBranchAddress("pz",&pz);
	const int events=tree->GetEntries();
	const int draw_events = 10;
	int smaller=min(events,draw_events);
	const int chists = smaller;
	TH2F *hevent[chists];
	for(int m=0;m<chists;m++) {
		hevent[m] = new TH2F(Form("hevent_%d",m),"",100,-100,100,250,-50,200);//xz
	}
	cout<<"Loop over "<<events<<" events and draw "<<smaller<<" tracks"<<endl;
	for(int i=0;i<events;i++) {
		tree->GetEntry(i);
// 		cout<<px->size()<<endl;
		for(int j=0;j<px->size();j++) {
			if(i<draw_events) hevent[i]->Fill(px->at(j),pz->at(j),pE->at(j));
//			cout<<px->at(j)<<"\t"<<pz->at(j)<<endl;
		}
		if(i<draw_events) hevent[i]->Draw("COLZ");
	}
	TFile *fout = new TFile("tracks.root","recreate");
	for(int m=0;m<chists;m++) hevent[m]->Write();
	fout->Close();
}
