void makehist() {
	ifstream in;
	in.open("out7p5MeV.txt");
	TH1F *h1 = new TH1F("h1","",400,4,10);
	TRandom3 *rndm3 = new TRandom3();
	while(in.good()) {
		double Edep;
		in>>Edep;
		double ESi=rndm3->Gaus(Edep,0.04);//40 keV sigma
		h1->Fill(ESi);
	}
	h1->Draw();
}
