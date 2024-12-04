TGraph *gXS_0 = new TGraph("XS_0.txt");
TGraph *gXS_1 = new TGraph("XS_1.txt");
TGraph *gXS_2 = new TGraph("XS_2.txt");
TGraph *gstopping = new TGraph("stopping_power.txt");
TGraph *gXS_t = new TGraph("XS_t.txt"); //Atomic Data and Nuclear Data Tables, Vol. 15, No. 1, January 1975
/*
TGraph2D *g_gs = new TGraph2D("n0_ang.dat");
TGraph2D *gp  = new TGraph2D("n1_ang.dat");
TGraph2D *gpp = new TGraph2D("n2_ang.dat");
//*/
/*
TString filename_n0 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n0_ang-iso.dat";
TString filename_n1 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n1_ang-iso.dat";
TString filename_n2 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n2_ang-iso.dat";
//*/
///*
TString filename_n0 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n0_ang.dat";
TString filename_n1 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n1_ang.dat";
TString filename_n2 = "/home/filippo/Documents/Geant4/Simulations/AmBe/AmBeData/LegendrePols/n2_ang.dat";
//*/
TGraph2DErrors *g_gs;
TGraph2DErrors *gp;
TGraph2DErrors *gpp;
/*
TGraph2D *g_gs = new TGraph2D("n0_ang-iso.dat");
TGraph2D *gp  = new TGraph2D("n1_ang-iso.dat");
TGraph2D *gpp = new TGraph2D("n2_ang-iso.dat");
//*/

#include <CLHEP/Random/RandFlat.h>
#define G4UniformRand() CLHEP::HepRandom::getTheEngine()->flat()

TRandom3 *rndm = new TRandom3();
double d2r=atan(1.)/45.;
double pi=4.*atan(1.);
double mb=4.002603;//Mass 4He
double mt=9.012182;//Mass 9Be
double ml=1.0008665;//Mass neutron
double mh=12.;//Mass 12C
double amu=931.4941;
double Q=+5.702;
double max_depth = 0.03639234;//mm
double InteractionE(double);
double InteractionTheta(double);
void makesumgraph();
double Eloss(double,double);
double ConvertCMtoLab(double Eint,double thetaCM);
double NeutronEnergy(double Eint, double thetalab, double thetaCM);
double Einitial=5.48556;//MeV
double Eth=0. ;//Threshold energy
double thetaCM,thetalab;
TH1F *h1 = new TH1F("h1",";Neutron Energy [MeV];Arbitrary Units/200 keV bin",75,0,15);
TH1F *hEn = new TH1F("hEn",";Neutron Energy [MeV];Yield",1500,0.0,15.0);
TH1F *hth = new TH1F("hth","",1800,0.0,180);
TH2F *hEth = new TH2F("hEth","",180,0,180,1500,0,15.);
TH2F *hEthCM = new TH2F("hEthCM","",180,0,180,1500,0,15.);
TH2F *hthth = new TH2F("hthth","COM vs lab",1800,0,180,1800,0,180);
bool secondks=false;
bool excited=false;
bool excited2=false;
double GetNeutron(double alphaEnergy);

/*
void ambe() {
	gStyle->SetOptStat(0);
	const int nsamples = 1000000;
	for(int n=0;n<nsamples;n++) {
		Einitial = 5.48556;//Main peak
		if(rndm->Rndm()<0.131) Einitial = 5.4428;//Satellite peak
		//Sample Ep_interaction
		double Eint = InteractionE(Einitial);
		h1->Fill(Eint);
		//Sample theta CM
		thetaCM = InteractionTheta(Eint);
		if(thetaCM==0) continue;//Skip this guy?
		thetalab = ConvertCMtoLab(Eint,thetaCM);
		hth->Fill((thetalab/d2r));
		double En = NeutronEnergy(Eint,thetalab,thetaCM);
		if(En==0) continue;//Skip this guy
		if(!excited2) continue;
		hthth->Fill(thetaCM/d2r,thetalab/d2r);
		hEn->Fill(En);
		hEth->Fill(thetalab/d2r,En);
//		if(!secondks)hEth->Fill(thetalab/d2r,En);
		hEthCM->Fill(thetaCM/d2r,En);
	}
	TCanvas *c = new TCanvas();
	hth->Draw();
	TCanvas *c2 = new TCanvas();
	hEn->Draw();
	TLegend *leg = new TLegend(0.7,0.9,0.7,0.9);
	leg->AddEntry(hEn,"BAmBe");
	leg->Draw("SAME");
	TCanvas *c3 = new TCanvas();
	c3->Divide(2,1);
	c3->cd(1);
	hEth->Draw("COLZ");
	c3->cd(2);
	hEthCM->Draw("COLZ");
	TCanvas *c4 = new TCanvas();
	hthth->Draw("COLZ");
	TCanvas *c5 = new TCanvas();
//	cout<<gXS_0->Eval(1.94999)<<endl;
	gXS_t->Draw();
	gXS_0->SetLineColor(kRed);
	gXS_0->Draw("SAME");
	gXS_1->SetLineColor(kBlue);
	gXS_1->Draw("SAME");
	gXS_2->SetLineColor(kGreen);
	gXS_2->Draw("SAME");
	TCanvas *c6 = new TCanvas();
	h1->Draw();
}
*/

void loadFile(TString filename, std::vector<double>& Z_array, std::vector<double>& X_array, std::vector<double>& Y_array, std::vector<double>& ex_array, std::vector<double>& ey_array)
{
  std::ifstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::string line;
  while (std::getline(file, line))
  {
    std::stringstream ss(line);
    double Z, X, Y, ex, ey;

    if (ss >> Z >> X >> Y >> ex >> ey)
    {
      Z_array.push_back(Z);
      X_array.push_back(X);
      Y_array.push_back(Y);
      ex_array.push_back(ex);
      ey_array.push_back(ey);
    }
    else
    {
      std::cerr << "Warning: Could not parse line: " << line << std::endl;
    }
  }

  file.close();
}

TGraph2DErrors* loadAng(TString filename)
{
    std::vector<double> Z_array, X_array, Y_array, ex_array, ey_array, ez_array;
    loadFile(filename, Z_array, X_array, Y_array, ex_array, ey_array);
	for (int i=0; i<Z_array.size(); i++) ez_array.push_back(0);

    // Display loaded data
	TGraph2DErrors *tg = new TGraph2DErrors(Z_array.size(),&Z_array[0],&X_array[0],&Y_array[0],&ez_array[0],&ex_array[0],&ey_array[0]);
	return tg;
}


void ambe()
{
  g_gs = loadAng(filename_n0);
  g_gs->SetName("n0_ang");
  gp = loadAng(filename_n1);
  gp->SetName("n1_ang");
  gpp = loadAng(filename_n2);
  gpp->SetName("n2_ang");

  double Einitial = 5.48556;//Main peak

  gStyle->SetOptStat(0);
  //const int nsamples = 1000000;
  const int nsamples = 2.27e6 * 1;
  for(int n=0;n<nsamples;n++)
  {
    //if(rndm->Rndm()<0.1338) Einitial = 5.44280;//Satellite peak
    if(G4UniformRand()<0.1338) Einitial = 5.44280;//Satellite peak
    double alphaEnergy = InteractionE(Einitial);
    h1->Fill(GetNeutron(alphaEnergy));
  }

  TCanvas *c1 = new TCanvas();
  h1->Draw();
}

double GetNeutron(double alphaEnergy)
{
  double Eint = alphaEnergy;//energy of interacting alpha
  //Sample theta CM
  thetaCM = InteractionTheta(Eint);
  if(thetaCM==0) return 0;//Skip this guy?
  thetalab = ConvertCMtoLab(Eint,thetaCM);
  double En = NeutronEnergy(Eint,thetalab,thetaCM);
  if(En==0) return 0;//Skip this guy
  return En;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double InteractionE(double Ep) {
	//Rejection sampling for z_interaction based on XS(Ep) and stopping power
	//Then return Ep for z
	bool trials=true;
	double XSmax = 700;
	int loopcounter=0;
	while(trials) {//Randomly choose a depth into the target (rather than energy to keep same target thickness)
    //double trial_z = max_depth*rndm->Rndm();
    double trial_z = max_depth*G4UniformRand();
		//Calc XS at Z
		double energy = Eloss(Ep,trial_z);
		if(energy<=Eth) continue;

    double XS = (gXS_0->Eval(energy)+gXS_1->Eval(energy)+gXS_2->Eval(energy));
    //cout<<gXS_1->Eval(energy)<<"\t"<<gXS_0->Eval(energy)<<endl;
    double choose_chan = G4UniformRand();
    double BRn0=gXS_0->Eval(energy)/XS;
    double BRn1=gXS_1->Eval(energy)/XS;
    double BRn2=gXS_2->Eval(energy)/XS;
		excited=false;
		excited2=false;
		if(choose_chan<BRn0) {
//			cout<<"Excited"<<endl;
			excited=false;
			excited2=false;
		}
		else if(choose_chan<BRn0+BRn1) {
			excited=true;
			excited2=false;
		}
		else {
			excited=false;
			excited2=true;
		}
//		cout<<energy<<"\t"<<XS<<endl;
		if(XS>XSmax) {cout<<"Max XS is exceeded!!"<<endl;}//
		//if(XS>XSmax*rndm->Rndm()) {
    if(XS>XSmax*G4UniformRand()) {
			return energy;
		}
		loopcounter++;
		if(loopcounter>1000) {
			cout<<"Rejection sampling not working properly - 1000 samples taken"<<endl;
			return 0;
		}
	}
	return 0;
}

/*
double Eloss(double Ea, double z) {//Energy of alpha of initial energy Ea after z mm in Be
	double range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
//	cout<<"INitial range\t"<<range<<endl;
	if(z>range) return 0;//Fully stopped before Z
	range-=z;
	double Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
//	cout<<Ep<<"\t"<<range<<"\t"<<z<<"\t"<<Eout<<endl;
	return Eout;
}
*/
double Eloss(double Ea, double z) //Energy of alpha of initial energy Ea after z mm in Be
{
  /*
   * Manually extracted from Lise++
   * Using 241Am at stechio 110
   * 16O at stechio 220
   * 9Be at stechio 890
   *
   * and alpha energy at max-> 5.48556
   *
   * -> max range 36.82486*um;
  */
  //double range=0.000467*Ea*Ea+0.00238*Ea+0.000481;//mm
  double range=0.000528638*Ea*Ea+0.00355149*Ea+0.000974915;//mm
  if(z>range) return 0;//Fully stopped before Z
  range-=z;
  //double Eout=-3232.85*range*range+282.754*range+0.0932;//MeV
  double Eout=-1775.23*range*range+216.585*range-0.0837652;//MeV
  return Eout;
}

double InteractionTheta(double Ep) //COM theta
{
  //return (100+10*G4UniformRand())*d2r;
  double random_theta = acos(-1+2*G4UniformRand());
  //random_theta = 2.*pi*G4UniformRand();//NOTE:check this
  random_theta = pi*G4UniformRand();
  double A[3]={1,0,0};
  bool sample=true;
  double XSmax=1.2;//Max is 1, XS data are in realtive yield
  int counter=0;
  TGraph2DErrors *select;
  if(!excited && !excited2) {select=g_gs;}
  if(excited && !excited2) {select=gp;}
  if(!excited && excited2) {select=gpp;}
  while(sample)
  {
    double XS_sample = XSmax*G4UniformRand();
    random_theta = acos(-1+2.*G4UniformRand());
    double XS=select->Interpolate(Ep,random_theta/d2r);
    if(XS>XSmax) std::cout<<"Sampling error for angle"<<std::endl;
    if(XS==0 || excited2) XS=1;//Outside of the energy range for this data set or Hoyle (poor data)
    if(XS_sample<XS)
    {
      return random_theta;//Isotropic for now
    }
    counter++;
    if(counter>100)
    {
      std::cout<<"Over 100 samples taken - no solution for Ealpha = "<<Ep<<"\t"<<excited<<"\t"<<excited2<<std::endl;
      return 0;
    }
  }
  return -1;//should not reach //TODO:check if correct
}

double ConvertCMtoLab(double Ep, double thetaCM) {//Convert thetaCM to lab
	double Ex=0;
	if(excited) Ex=4.44;
	if(excited2) Ex=7.654;

	double gamma = sqrt(((mb*ml)/(mt*mh))*(Ep/(Ep+(Q-Ex)*(1.+mb/mt))));//Calculate gamma for the conversion to lab COM
	gamma = gamma*((mb+mt)/(ml+mh));//Account for change in the CM velocity
	double thetalab = atan2(sin(thetaCM),(cos(thetaCM)+gamma));//
	secondks = false;
	if(gamma*cos(thetaCM)<=-1.) secondks=true;//Neutron is going bac
//	thetalab=(12*rndm->Rndm())*d2r;//TEST
//	if(rndm->Rndm()>0.5) secondks=true;//TEST
//	cout<<secondks<<"\t"<<thetaCM/d2r<<"\t"<<gamma*cos(thetaCM)<<"\t"<<gamma<<endl;
	return thetalab;
}

double NeutronEnergy(double Eb, double theta, double thetaCM) {//Get neutron energy for a given Ep and theta
	double En;
	double Ex=0;
	if(excited) Ex=4.44;//MeV
	if(excited2) Ex=7.654;

	double ECM = mt*Eb/(mt+mb)+Q-Ex;
//	cout<<ECM<<endl;
	if(ECM<0) {
		cout<<"Below threshold! Something is wrong!"<<endl;
		return 0;
	}
	double p_nT = sqrt(2.*ECM*ml*(mh/(ml+mh)));
	double phi=2*pi*G4UniformRand();
	double pn[3]={p_nT*sin(thetaCM)*cos(phi),p_nT*sin(thetaCM)*sin(phi),p_nT*cos(thetaCM)};//COM neutron momentum
	pn[2]+=(ml/(mh+ml))*sqrt(2.*mb*Eb);//Add the momentum boost for the lab frame
	double En_new = 0;
	for(int i=0;i<3;i++) En_new+=pn[i]*pn[i];
	En_new /= (2.*ml);
	En = En_new;
	theta = acos(pn[2]/p_nT);
	return En;
}
