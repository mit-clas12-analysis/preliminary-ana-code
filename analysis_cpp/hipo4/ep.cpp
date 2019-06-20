#include <ep.h>
// Forward-declaring functions

// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV (proton mass      )
	double mPiC    = 0.13957; //GeV (charged pion mass)
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;

	// ----------------------------------------------------------------------------------
	// Getting input arguments
	TString inputFile;
	double Ebeam, mtar;

	if(argc==3){
		if(atoi(argv[1])==1){
			cout << "Will assume this hipo file corresponds to: Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
			Ebeam = 6.4; //GeV
			mtar  = mp;
		}
		else if(atoi(argv[1])==2){
			cout << "Will assume this hipo file corresponds to: Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
			Ebeam = 10.6; //GeV
			mtar  = mD;
		}
		inputFile = argv[2];
	}
	else {
		cout << "=========================\nRun this code as:\n./code A path/to/input/file\n" << endl;
		cout << "where: A = 1 -> Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
		cout << "         = 2 -> Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
		cout << "=========================" << endl;
		exit(0);
	}

	TVector3 V3_Ebeam(0,0,Ebeam);
	TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);
	TLorentzVector V4_mtar(0,0,0,mtar);

	// ----------------------------------------------------------------------------------
	// Event selection cuts
	double cut_ep      =     2; //GeV
	double cut_chi2pid =     5;
	double cut_min_vz  =   -15; //cm
	double cut_max_vz  =    10; //cm
	double cut_W       =     0; //GeV
	double cut_uvw     =    15; //cm
	double cut_Epcal   = 0.060; //GeV (60 MeV)
	double cut_tof_e   =    10; //ns
	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
   hipo::reader  reader;
   reader.open(inputFile);

   hipo::dictionary  factory;

   reader.readDictionary(factory);

   factory.show();
   hipo::structure  particles;
   hipo::structure  detectors;

   hipo::event      event;
   int counter = 0;

   hipo::bank  dataPART;
   hipo::bank  dataCALO;

   hipo::bank PART(factory.getSchema("REC::Particle"));
  

   while(reader.next()==true){
      reader.read(event);
      //reader.next();
      //event.getStructure(dataBank,30,1);
      //dataBank.show();
      event.getStructure(PART);
      //PART.show();
      int nrows = PART.getRows();

      for(int i = 0; i < nrows; i++){
        int   pid = PART.getInt("pid",i);
        float  px = PART.getFloat("px",i);
        float  py = PART.getFloat("py",i);
        float  pz = PART.getFloat("pz",i);

        /*int   pid = PART.getInt(1,i);
        float  px = PART.getFloat(2,i);
        float  py = PART.getFloat(3,i);
        float  pz = PART.getFloat(4,i);*/
        printf("%6d %8.4f %8.4f %8.4f\n",pid,px,py,pz);
      }
      //event.getStructure(dataCALO,300,32);
      //dataPART.show();
      counter++;
   }
   printf("processed events = %d\n",counter);

    // // processEvent(event);
	// // particle     particles   ("REC::Particle"    ,reader);
	// // clas12calorimeter  calo        ("REC::Calorimeter" ,reader);
	// // BScintillator scintillator("REC::Scintillator",reader);

	// int event_counter = 0;
	// // ----------------------------------------------------------------------------------
	// // Loop over events and print them on the screen
	// while(reader.next()==true){

	// 	if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
	// 	event_counter++;

	// 	//particles.show();
	// 	//calo.show();
	// 	//event.show();
	// 	//bank_scintillator.show();	

	// 	// Particle bank
	// 	int pid0       = particles.getPid    (0);	// electron candidate id assigned by clas
	// 	TVector3 V3_ev = particles.getV3v    (0);	// electron candidate vertex vector
	// 	TVector3 V3_ep = particles.getV3P    (0);	// electron candidate momentum vector
	// 	float chr0     = particles.getCharge (0);	// electron candidate charge
	// 	float eBeta    = particles.getBeta   (0);	// electron candidate beta = v/c
	// 	float chi2pid  = particles.getChi2pid(0);	// electron candidate goodness of pid fit
	// 	int eStatus    = particles.getStatus (0);	// electron candidate status

	// 	// Calorimeter bank	
	// 	float Epcal = calo.getPcalE(0); 
	// 	float Ee    = calo.getTotE (0);
	// 	float lU    = calo.getLU   (0);	// electron candidate distance on U-side [cm?]
	// 	float lV    = calo.getLV   (0);	// electron candidate distance on V-side [cm?]
	// 	float lW    = calo.getLW   (0);	// electron candidate distance on W-side [cm?]

	// 	if(Ee==0) continue;

	// 	// Event bank
	// 	double t_vtx   = event.getSTTime(0);

	// 	// Scintillator bank
	// 	double t_e     = scintillator.getTime(0);

	// 	// calculated variables
	// 	double ep     = V3_ep.Mag();		// electron candidate momentum magnitude [GeV]
	// 	TLorentzVector V4_ep(V3_ep,ep);
	// 	double tof_e  = t_e - t_vtx;		// electron candidate time-of-flight [ns]

	// 	// Transfer variables
	// 	TVector3 V3_q = V3_Ebeam - V3_ep;
	// 	double Q2     = 4*ep*Ebeam*pow(TMath::Sin(V3_ep.Theta()/2.),2);       // Q-squared [GeV^2]
	// 	double omega  = Ebeam - ep;                                              // Transfer energy [GeV]
	// 	double W2     = mp*mp-Q2+2*omega*mp;
	// 	double xB     = Q2/2./mp/omega; 

	// 	TLorentzVector V4_q(V3_q,omega);

	// 	// -------------------------------------------------------------------------
	// 	// Fill some histograms before cutting on good electrons
	// 	h2_e_Ep_p_0 -> Fill(ep,Ee/ep);
	// 	h2_e_vz_phi -> Fill(rad2deg*V3_ep.Phi(),V3_ev.Z());

	// 	// -------------------------------------------------------------------------
	// 	// Electron PID from Dan Carman
	// 	// - (DONE)	pid=11 from EB
	// 	// - (DONE)	p > 2 GeV
	// 	// - (DONE)	p < Ebeam
	// 	// - (DONE)	TOF > 10 ns //(need bank 330) event.json
	// 	// - (DONE)	vz: -15 to 10 cm
	// 	// - (DONE)	W > 0 GeV
	// 	// - 		Sampling fraction +/-3sigma about E/p vs. p mean
	// 	// - (DONE)	~15 cm fiducial cuts on U, V, W to contain full shower (need bank 332 lu, lv, lw)
	// 	// - (DONE)	abs(chisq PID) < 5 (goodness of PID from EB)
	// 	// - (DONE)	PCAL > 60 MeV (to remove min-i) (bank 332 layers 1(PCAL), 4(EC inner), 7(EC outter))
	// 	// -------------------------------------------------------------------------
	// 	// Only keep events for which the first particle is an electron
	// 	if(             (pid0!=11              )||
	// 			(chr0!=-1              )||
	// 			(chi2pid>=cut_chi2pid  )||
	// 			(ep<=cut_ep            )||
	// 			(ep>=Ebeam             )||
	// 			(V3_ev.Z()>cut_max_vz  )||
	// 			(V3_ev.Z()<cut_min_vz  )||
	// 			(lU<cut_uvw            )||
	// 			(lV<cut_uvw            )||
	// 			(lW<cut_uvw            )||
	// 			(Epcal<cut_Epcal       )||
	// 			(TMath::Sqrt(W2)<=cut_W)||
	// 			(tof_e<cut_tof_e       )
	// 	  ) continue;

	// 	// -------------------------------------------------------------------------
	// 	// Filling electron histograms
	// 	h1_e_lu  -> Fill(lU                   );
	// 	h1_e_lv  -> Fill(lV                   );
	// 	h1_e_lw  -> Fill(lW                   );
	// 	h1_e_px  -> Fill(V3_ep.X()            );
	// 	h1_e_py  -> Fill(V3_ep.Y()            );
	// 	h1_e_pz  -> Fill(V3_ep.Z()            );
	// 	h1_e_p   -> Fill(V3_ep.Mag()          );
	// 	h1_e_vz  -> Fill(V3_ev.Z()            );
	// 	h1_W     -> Fill(TMath::Sqrt(W2)      );
	// 	h1_xB    -> Fill(xB                   );
	// 	h1_e_th  -> Fill(rad2deg*V3_ep.Theta());
	// 	h1_e_phi -> Fill(rad2deg*V3_ep.Phi()  );
	// 	h1_e_tof -> Fill(tof_e                );

	// 	h2_e_Ep_p_1 -> Fill(ep           ,Ee/ep       );
	// 	h2_e_th_phi -> Fill(rad2deg*V3_ep.Phi(),rad2deg*V3_ep.Theta());
	// 	h2_e_tof_p  -> Fill(ep           ,tof_e       );
	// 	h2_pe_the   -> Fill(rad2deg*V3_ep.Theta(),ep  );

	// 	if((xB>1.2)||(xB<0.8)) continue;

	// 	// -------------------------------------------------------------------------

	// 	// pi+ variables
	// 	int nProtons = 0;
	// 	int tmp_fast_p_idx  = 0;	// index of the pi+

	// 	int nParticles = particles.getSize();	

	// 	for(int par = 1; par < nParticles; par++) {
	// 		// Particle bank

	// 		int pidi       = particles.getPid    (par);       // id assigned by clas
	// 		TVector3 V3_hv = particles.getV3v    (par);       // vertex vector
	// 		TVector3 V3_hp = particles.getV3P    (par);       // momentum vector
	// 		float chri     = particles.getCharge (par);       // charge
	// 		float beta_p   = particles.getBeta   (par);       // beta = v/c
	// 		double pp     = V3_hp.Mag();

	// 		if     (chri== 1) h2_beta_p_pos -> Fill(V3_hp.Mag(), beta_p  );

	// 		// -------------------------------------------------------------------------
	// 		// proton
	// 		if(             (pidi==2212         )&&
	// 				(chri== 1           )&&
	// 				(V3_hp.Mag()<Ebeam  )
	// 		  ) {
	// 			nProtons++;
	// 			tmp_fast_p_idx=par;
	// 		}

	// 		// ----------------------------------------------------------------------

	// 	}

	// 	if(nProtons==1&&nParticles==2) {

	// 		TVector3 V3_pv  = particles.getV3v  (tmp_fast_p_idx);		// p 2 vertex vector [cm]
	// 		TVector3 V3_pp  = particles.getV3P  (tmp_fast_p_idx);		// p 2 momentum vector [GeV]
	// 		float beta_p    = particles.getBeta (tmp_fast_p_idx);		// p 2 beta = v/c

	// 		TLorentzVector V4_pp;
	// 		V4_pp.SetXYZM( V3_pp.X() , V3_pp.Y() , V3_pp.Z() , mp );

	// 		// Missing momentum components
	// 		TVector3       V3_Pm = V3_pp - V3_q;
	// 		TLorentzVector V4_Pm = V4_Ebeam + V4_mtar - V4_ep - V4_pp;	

	// 		double Mmiss = V4_Pm.M();		

	// 		double Em = fn_Emiss( V3_Pm.Mag() , omega , mtar , V4_pp.E() , mp );

	// 		// Predicted proton momentum and angle
	// 		double th_p_calc = TMath::ACos((Ebeam - ep*TMath::Cos(V3_ep.Theta()))/V3_pp.Mag());
	// 		double p_p_calc  = TMath::Sqrt(pow(Ebeam+mp-ep, 2)-mp*mp);

	// 		if(V3_Pm.Mag()>0.25&&Em<0.1){
	// 			// Filling histograms
	// 			h1_p_vz        -> Fill(V3_pp.Z()              );
	// 			h1_dlt_vz_ep   -> Fill(V3_pp.Z()  - V3_ev.Z() );
	// 			h1_p_px        -> Fill(V3_pp.X()              );
	// 			h1_p_py        -> Fill(V3_pp.Y()              );
	// 			h1_p_pz        -> Fill(V3_pp.Z()              );
	// 			h1_p_p         -> Fill(V3_pp.Mag()            );
	// 			h1_p_th        -> Fill(rad2deg*V3_pp.Theta()  );
	// 			h1_p_phi       -> Fill(rad2deg*V3_pp.Phi()    );
	// 			h1_pmx         -> Fill(V3_Pm.X()                );
	// 			h1_pmy         -> Fill(V3_Pm.Y()                );
	// 			h1_pmz         -> Fill(V3_Pm.Z()                );
	// 			h1_pm_th       -> Fill(rad2deg*V3_Pm.Theta()    );
	// 			h1_pm_ph       -> Fill(rad2deg*V3_Pm.Phi()      );
	// 			h1_pmiss       -> Fill(V3_Pm.Mag()              );
	// 			h1_Mmiss       -> Fill(Mmiss                    );

	// 			h1_p_th_meas_calc -> Fill((V3_pp.Theta()-th_p_calc)*rad2deg);
	// 			h1_p_p_meas_calc  -> Fill(V3_pp.Mag()-p_p_calc    );

	// 			h2_p_th_phi   -> Fill(rad2deg*V3_pp.Phi(), rad2deg*V3_pp.Theta());
	// 			h2_p_vz_phi   -> Fill(rad2deg*V3_pp.Phi(), V3_pp.Z()            );
	// 			h2_beta_p_p   -> Fill(V3_pp.Mag()        , beta_p               );
	// 			h2_pe_pp      -> Fill(V3_pp.Mag()        , ep                     );
	// 		}
	// 	}


	// 	h1_p_num -> Fill((double)(nProtons));

	// } // End of while loop over events

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2,TString titx,TString tity) {
	h2 -> GetXaxis() -> SetTitle(titx);
	h2 -> GetYaxis() -> SetTitle(tity);
}
// ========================================================================================================================================
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - TMath::Sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;                                                             			// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;
}

void InitiateHistograms(){
	// ----------------------------------------------------------------------------------
	// Declaring histograms
	// 1D histograms
	h1_e_vz  = new TH1F("h1_e_vz"  ,"h1_e_vz"  ,100, -50, 50);	PrettyTH1F(h1_e_vz  ,"v_{z} [cm]"           ,"Counts",62);
	h1_e_tof = new TH1F("h1_e_tof" ,"h1_e_tof" ,100,  20, 27);	PrettyTH1F(h1_e_tof ,"electron TOF [ns]"    ,"Counts",62);
	h1_e_px  = new TH1F("h1_e_px"  ,"h1_e_px"  ,100,  -2,  2);	PrettyTH1F(h1_e_px  ,"electron p_{x} [GeV]" ,"Counts",62);
	h1_e_py  = new TH1F("h1_e_py"  ,"h1_e_py"  ,100,  -2,  2);	PrettyTH1F(h1_e_py  ,"electron p_{y} [GeV]" ,"Counts",62);
	h1_e_pz  = new TH1F("h1_e_pz"  ,"h1_e_pz"  ,100,   0, 10);	PrettyTH1F(h1_e_pz  ,"electron p_{z} [GeV]" ,"Counts",62);
	h1_e_p   = new TH1F("h1_e_p"   ,"h1_e_p"   ,100,   0, 10);	PrettyTH1F(h1_e_p   ,"electron |p| [GeV]"   ,"Counts",62);
	h1_e_th  = new TH1F("h1_e_th"  ,"h1_e_th"  ,100,   0, 30);	PrettyTH1F(h1_e_th  ,"#theta_e [deg]"       ,"Counts",62);
	h1_e_phi = new TH1F("h1_e_phi" ,"h1_e_phi" ,100,-190,190);	PrettyTH1F(h1_e_phi ,"#phi_e [deg]"         ,"Counts",62);
	h1_e_lu  = new TH1F("h1_e_lu"  ,"h1_e_lu"  ,100,   0,500);	PrettyTH1F(h1_e_lu  ,"distance on U-side"   ,"Counts",62);
	h1_e_lv  = new TH1F("h1_e_lv"  ,"h1_e_lv"  ,100,   0,500);	PrettyTH1F(h1_e_lv  ,"distance on V-side"   ,"Counts",62);
	h1_e_lw  = new TH1F("h1_e_lw"  ,"h1_e_lw"  ,100,   0,500);	PrettyTH1F(h1_e_lw  ,"distance on W-side"   ,"Counts",62);
	h1_p_vz  = new TH1F("h1_p_vz"  ,"h1_p_vz"  ,100, -50, 50);	PrettyTH1F(h1_p_vz  ,"v_{z} [cm]"           ,"Counts",2);
	h1_p_num = new TH1F("h1_p_num" ,"h1_p_num" , 20,   0, 10);	PrettyTH1F(h1_p_num ,"#pi+ number"          ,"Counts",62);
	h1_p_px  = new TH1F("h1_p_px"  ,"h1_p_px"  ,100,  -2,  2);	PrettyTH1F(h1_p_px  ,"#pi+ p_{x} [GeV]"     ,"Counts",62);
	h1_p_py  = new TH1F("h1_p_py"  ,"h1_p_py"  ,100,  -2,  2);	PrettyTH1F(h1_p_py  ,"#pi+ p_{y} [GeV]"     ,"Counts",62);
	h1_p_pz  = new TH1F("h1_p_pz"  ,"h1_p_pz"  ,100,   0,  9);	PrettyTH1F(h1_p_pz  ,"#pi+ p_{z} [GeV]"     ,"Counts",62);
	h1_p_p   = new TH1F("h1_p_p"   ,"h1_p_p"   ,100,   0,  9);	PrettyTH1F(h1_p_p   ,"#pi+ |p| [GeV]"       ,"Counts",62);
	h1_p_th  = new TH1F("h1_p_th"  ,"h1_p_th"  ,100,   0, 80);	PrettyTH1F(h1_p_th  ,"#theta_p [deg]"       ,"Counts",62);
	h1_p_phi = new TH1F("h1_p_phi" ,"h1_p_phi" ,100,-190,190);	PrettyTH1F(h1_p_phi ,"#phi_p [deg]"         ,"Counts",62);
	h1_pmiss = new TH1F("h1_pmiss" ,"P_{miss}" ,100,   0,  5);	PrettyTH1F(h1_pmiss ,"Pm [GeV]"             ,"Counts",62);
	h1_pmx   = new TH1F("h1_pmx"   ,"h1_pmx"   ,100,  -2,  2);	PrettyTH1F(h1_pmx   ,"Pmx [GeV]"            ,"Counts",62);
	h1_pmy   = new TH1F("h1_pmy"   ,"h1_pmy"   ,100,  -2,  2);	PrettyTH1F(h1_pmy   ,"Pmy [GeV]"            ,"Counts",62);
	h1_pmz   = new TH1F("h1_pmz"   ,"h1_pmz"   ,100,  -3,  3);	PrettyTH1F(h1_pmz   ,"Pmz [GeV]"            ,"Counts",62);
	h1_pm_th = new TH1F("h1_pm_th" ,"h1_pm_th" ,100,   0,180);	PrettyTH1F(h1_pm_th ,"#theta_{Pm} [deg]"    ,"Counts",62);
	h1_pm_ph = new TH1F("h1_pm_ph" ,"h1_pm_ph" ,100,-190,190);	PrettyTH1F(h1_pm_ph ,"#phi_{Pm} [deg]"      ,"Counts",62);
	h1_Mmiss = new TH1F("h1_Mmiss" ,"h1_Mmiss" ,100,   0,  4);	PrettyTH1F(h1_Mmiss ,"m_{miss} [GeV]"       ,"Counts",62);
	h1_Em    = new TH1F("h1_Em"    ,"h1_Em"    ,100,  -3,  3);	PrettyTH1F(h1_Em    ,"Em [GeV]"             ,"Counts",62);
	h1_W     = new TH1F("h1_W"     ,"h1_W"     ,100,   0,  4);	PrettyTH1F(h1_W     ,"W [GeV]"              ,"Counts",62);
	h1_xB    = new TH1F("h1_xB"    ,"h1_xB"    ,100,   0,  4);	PrettyTH1F(h1_xB    ,"x_B"                  ,"Counts",62);
	h1_dlt_vz_ep  = new TH1F("h1_dlt_vz_ep"  ,"electron - #pi+",100, -20, 20);	PrettyTH1F(h1_dlt_vz_ep  ,"#Delta v_{z} [cm]"    ,"Counts",2);
	h1_p_th_meas_calc = new TH1F("h1_p_th_meas_calc","h1_p_th_meas_calc",80,-50,50);
	h1_p_p_meas_calc = new TH1F("h1_p_p_meas_calc","h1_p_p_meas_calc",80,-5,5);

	PrettyTH1F(h1_p_th_meas_calc,"proton #theta_{measured}-#theta_{calculated} [deg]","Counts",62);
	PrettyTH1F(h1_p_p_meas_calc,"proton p_{measured}-p_{calculated} [GeV]"    ,"Counts",62);

	// 2D histograms
	h2_e_Ep_p_0   = new TH2F("h2_e_Ep_p_0"   ,"h2_e_Ep_p_0"  ,100,1   , 10,100,  0,0.4);
	h2_e_Ep_p_1   = new TH2F("h2_e_Ep_p_1"   ,"h2_e_Ep_p_1"  ,100,1   , 10,100,  0,0.4);	
	h2_e_th_phi   = new TH2F("h2_e_th_phi"   ,"h2_e_th_phi"  ,100,-190,190,100,  0, 30);
	h2_p_th_phi   = new TH2F("h2_p_th_phi"   ,"h2_p_th_phi"  ,100,-190,190,100,  0, 80);
	h2_beta_p_pos = new TH2F("h2_beta_p_pos" ,"h2_beta_p_pos",100,0   ,  6,100,0.1,1.1);
	h2_beta_p_p   = new TH2F("h2_beta_p_p"   ,"h2_beta_p_p"  ,100,0   ,  6,100,0.1,1.1);
	h2_e_vz_phi   = new TH2F("h2_e_vz_phi"   ,"h2_e_vz_phi"  ,100,-190,190,100,-50, 50);
	h2_p_vz_phi   = new TH2F("h2_p_vz_phi"   ,"h2_p_vz_phi"  ,100,-190,190,100,-50, 50);
	h2_e_tof_p    = new TH2F("h2_e_tof_p"    ,"h2_e_tof_p"   ,100,1   ,  7,100, 20, 27);
	h2_p_dtT_p_0  = new TH2F("h2_p_dtT_p_0"  ,"h2_p_dtT_p_0" ,100,0   ,  6,100,- 4,  4);
	h2_p_dtT_p_1  = new TH2F("h2_p_dtT_p_1"  ,"h2_p_dtT_p_1" ,100,0   ,  6,100,- 4,  4);
	h2_p_tof_det  = new TH2F("h2_p_tof_det"  ,"h2_p_tof_det" , 72,0   , 36,100,  0, 40);
	h2_p_dtT_det  = new TH2F("h2_p_dtT_det"  ,"h2_p_dtT_det" , 72,0   , 36,100,-20, 20);
	h2_Em_Pm      = new TH2F("h2_Em_Pm"      ,"h2_Em_Pm"     ,100,0   ,  3,100,-2 ,  2);
	h2_pe_pp      = new TH2F("h2_pe_pp"      ,"h2_pe_pp"     ,100,0   ,  6,100, 1 , 10);
	h2_pe_the     = new TH2F("h2_pe_the"     ,"h2_pe_the"    ,100,0   , 30,100, 1 , 10);

	PrettyTH2F(h2_e_Ep_p_0  ,"p_{e} [GeV]"   ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_Ep_p_1  ,"p_{e} [GeV]"   ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_th_phi  ,"#phi_e [deg]"  ,"#theta_e [deg]"      );
	PrettyTH2F(h2_p_th_phi  ,"#phi_p [deg]"  ,"#theta_p [deg]"      );
	PrettyTH2F(h2_beta_p_pos,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_beta_p_p  ,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_e_vz_phi  ,"#phi_e [deg]"  ,"e v_{z} [cm]"        );
	PrettyTH2F(h2_p_vz_phi  ,"#phi_p [deg]"  ,"p v_{z} [cm]"        );
	PrettyTH2F(h2_e_tof_p   ,"p_{e} [GeV]"   ,"electron TOF [ns]"   );
	PrettyTH2F(h2_p_dtT_p_0 ,"p_{e} [GeV]"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_dtT_p_1 ,"p_{e} [GeV]"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_tof_det ,"detector id"   ,"candidate p tof [ns]");
	PrettyTH2F(h2_p_dtT_det ,"detector id"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_Em_Pm     ,"Pm [GeV]"      ,"Em [GeV]"            );
	PrettyTH2F(h2_pe_pp     ,"p p [GeV]"     ,"e p [GeV]"           );
	PrettyTH2F(h2_pe_the    ,"#theta_e [deg]","p_{e} [GeV]"         );

		H_elec_theta_mom = new H2F("H_elec_theta_mom","H_elec_theta_mom",100,0,EB,100,0,40);
		H_elec_theta_mom.setTitle("elec #theta vs p");
		H_elec_theta_mom.setTitleX("p (GeV)");
		H_elec_theta_mom.setTitleY("#theta (^o)");
		H_elec_theta_phi = new H2F("H_elec_theta_phi","H_elec_theta_phi",100,-180,180,100,0,40);
		H_elec_theta_phi.setTitle("elec #theta vs #phi");
		H_elec_theta_phi.setTitleX("#phi (^o)");
		H_elec_theta_phi.setTitleY("#theta (^o)");
		H_elec_phi_mom = new H2F("H_elec_phi_mom","H_elec_phi_mom",100,0,EB,100,-180,180);
		H_elec_phi_mom.setTitle("elec #phi vs p");
		H_elec_phi_mom.setTitleX("p (GeV)");
		H_elec_phi_mom.setTitleY("#phi (^o)");
		H_prot_theta_mom = new H2F("H_prot_theta_mom","H_prot_theta_mom",100,0,11,100,0,90);
		H_prot_theta_mom.setTitle("prot #theta vs p");
		H_prot_theta_mom.setTitleX("p (GeV)");
		H_prot_theta_mom.setTitleY("#theta (^o)");
		H_prot_theta_phi = new H2F("H_prot_theta_phi","H_prot_theta_phi",100,-180,180,100,0,90);
		H_prot_theta_phi.setTitle("prot #theta vs #phi");
		H_prot_theta_phi.setTitleX("#phi (^o)");
		H_prot_theta_phi.setTitleY("#theta (^o)");
		H_prot_phi_mom = new H2F("H_prot_phi_mom","H_prot_phi_mom",100,0,EB,100,-180,180);
		H_prot_phi_mom.setTitle("prot #phi vs p");
		H_prot_phi_mom.setTitleX("p (GeV)");
		H_prot_phi_mom.setTitleY("#phi (^o)");
		H_g_theta_E = new H2F("H_g_theta_E","H_g_theta_E",100,0,EB,100,0,40);
		H_g_theta_E.setTitle("#gamma #theta vs p");
		H_g_theta_E.setTitleX("p (GeV)");
		H_g_theta_E.setTitleY("#theta (^o)");
		H_g_theta_phi = new H2F("H_g_theta_phi","H_g_theta_phi",100,-180,180,100,0,40);
		H_g_theta_phi.setTitle("#gamma #theta vs #phi");
		H_g_theta_phi.setTitleX("#phi (^o)");
		H_g_theta_phi.setTitleY("#theta (^o)");
		H_g_phi_E = new H2F("H_g_phi_E","H_g_phi_E",100,0,EB,100,-180,180);
		H_g_phi_E.setTitle("#gamma #phi vs p");
		H_g_phi_E.setTitleX("p (GeV)");
		H_g_phi_E.setTitleY("#phi (^o)");
		H_g1_theta_E = new H2F("H_g1_theta_E","H_g1_theta_E",100,0,EB,100,0,40);
		H_g1_theta_E.setTitle("#gamma1 #theta vs p");
		H_g1_theta_E.setTitleX("p (GeV)");
		H_g1_theta_E.setTitleY("#theta (^o)");
		H_g1_theta_phi = new H2F("H_g1_theta_phi","H_g1_theta_phi",100,-180,180,100,0,40);
		H_g1_theta_phi.setTitle("#gamma1 #theta vs #phi");
		H_g1_theta_phi.setTitleX("#phi (^o)");
		H_g1_theta_phi.setTitleY("#theta (^o)");
		H_g1_phi_E = new H2F("H_g1_phi_E","H_g1_phi_E",100,0,EB,100,-180,180);
		H_g1_phi_E.setTitle("#gamma1 #phi vs p");
		H_g1_phi_E.setTitleX("p (GeV)");
		H_g1_phi_E.setTitleY("#phi (^o)");
		H_g2_theta_E = new H2F("H_g2_theta_E","H_g2_theta_E",100,0,EB,100,0,40);
		H_g2_theta_E.setTitle("#gamma2 #theta vs p");
		H_g2_theta_E.setTitleX("p (GeV)");
		H_g2_theta_E.setTitleY("#theta (^o)");
		H_g2_theta_phi = new H2F("H_g2_theta_phi","H_g2_theta_phi",100,-180,180,100,0,40);
		H_g2_theta_phi.setTitle("#gamma2 #theta vs #phi");
		H_g2_theta_phi.setTitleX("#phi (^o)");
		H_g2_theta_phi.setTitleY("#theta (^o)");
		H_g2_phi_E = new H2F("H_g2_phi_E","H_g2_phi_E",100,0,EB,100,-180,180);
		H_g2_phi_E.setTitle("#gamma2 #phi vs p");
		H_g2_phi_E.setTitleX("p (GeV)");
		H_g2_phi_E.setTitleY("#phi (^o)");
		H_prot_beta_p = new H2F("H_prot_beta_p","H_prot_beta_p",100,0,4,100,0,1.2);
		H_prot_beta_p.setTitle("prot #beta vs p");
		H_prot_beta_p.setTitleX("p (GeV)");
		H_prot_beta_p.setTitleY("#beta");
		H_phot_b_p = new H2F("H_phot_b_p","H_phot_b_p",100,0.5,EB,100,0.75,1.25);
		H_phot_b_p.setTitle("#gamma #beta vs p");
		H_phot_b_p.setTitleX("p (GeV)");
		H_phot_b_p.setTitleY("#beta");
		H_prot_e_vz = new H2F("H_prot_e_vz","H_prot_e_vz",100,-25,25,100,-25,25);
		H_prot_e_vz.setTitle("prot vz vs elec vz");
		H_prot_e_vz.setTitleX("elec vz (cm)");
		H_prot_e_vz.setTitleY("prot vz (cm)");

		H_prot_miss_phi = new H2F("H_prot_miss_phi","H_prot_miss_phi",100,-180,180,100,-180,180);
		H_prot_miss_phi.setTitle("DVCS proton #phi vs miss #phi");
		H_prot_miss_phi.setTitleX("proton #phi (^o)");
		H_prot_miss_phi.setTitleY("miss #phi (^o)");
		H_phot_elec_angle = new H2F("H_phot_elec_angle","H_phot_elec_angle",100,-180,180,100,0,40);
	       	H_phot_elec_angle.setTitle("angle #gamma e vs #phi");
		H_phot_elec_angle.setTitleX("#phi Trento (^o)");
		H_phot_elec_angle.setTitleY("angle #gamma e (^o)");
		H_phot_virt_angle = new H2F("H_phot_virt_angle","H_phot_virt_angle",100,0,5,100,0,35);
		H_phot_virt_angle.setTitle("angle #gamma#gamma^* vs -t");
		H_phot_virt_angle.setTitleX("-t (GeV^2)");
		H_phot_virt_angle.setTitleY("angle #gamma #gamma^* (^o)");

		H_prot_pi0_miss_phi = new H2F("H_prot_pi0_miss_phi","H_prot_pi0_miss_phi",100,-180,180,100,-180,180);
		H_prot_pi0_miss_phi.setTitle("#pi^0 proton #phi vs miss #phi");
		H_prot_pi0_miss_phi.setTitleX("proton #phi (^o)");
		H_prot_pi0_miss_phi.setTitleY("miss #phi (^o)");
		H_pi0_elec_angle = new H2F("H_pi0_elec_angle","H_pi0_elec_angle",100,-180,180,100,0,40);
	       	H_pi0_elec_angle.setTitle("angle #pi^0 e vs #phi");
		H_pi0_elec_angle.setTitleX("#phi Trento (^o)");
		H_pi0_elec_angle.setTitleY("angle #pi^0 e (^o)");
		H_pi0_virt_angle = new H2F("H_pi0_virt_angle","H_pi0_virt_angle",100,0,5,100,0,35);
		H_pi0_virt_angle.setTitle("angle #pi^0#gamma^* vs -t");
		H_pi0_virt_angle.setTitleX("-t (GeV^2)");
		H_pi0_virt_angle.setTitleY("angle #pi^0 #gamma^* (^o)");

		H_dvcs_Q2xB = new H2F("H_dvcs_Q2xB","H_dvcs_Q2xB",100,0,1,100,0,12);
		H_dvcs_Q2xB.setTitle("DVCS Q2 vs xB");
		H_dvcs_Q2xB.setTitleX("xB");
		H_dvcs_Q2xB.setTitleY("Q2");
		H_dvcs_tphi = new H2F("H_dvcs_tphi","H_dvcs_tphi",100,-180,180,100,0,5);
		H_dvcs_tphi.setTitle("DVCS -t vs #phi");
		H_dvcs_tphi.setTitleX("#phi (^o)");
		H_dvcs_tphi.setTitleY("-t (GeV^2)");
		H_dvcs_pT = new H1F("H_dvcs_pT","H_dvcs_pT",100,0,1.0);
		H_dvcs_pT.setTitle("DVCS missing pT");
		H_dvcs_pT.setTitleX("pT (GeV)");
		H_dvcs_cone = new H1F("H_dvcs_cone","H_dvcs_cone",100,0,5);
		H_dvcs_cone.setTitle("DVCS cone angle");
		H_dvcs_cone.setTitleX("cone ang (^o)");
		H_dvcs_copl = new H1F("H_dvcs_copl","H_dvcs_copl",100,0,30);
		H_dvcs_copl.setTitle("DVCS copl. angle");
		H_dvcs_copl.setTitleX("copl angle (^o)");
		H_dvcs_ME = new H1F("H_dvcs_ME","H_dvcs_ME",100,-1,2);
		H_dvcs_ME.setTitle("DVCS miss E");
		H_dvcs_ME.setTitleX("miss E (GeV)");
		H_dvcs_MM_ep = new H1F("H_dvcs_MM_ep","H_dvcs_MM_ep",100,-1,1.2);
		H_dvcs_MM_ep.setTitle("DVCS ep MM^2");
		H_dvcs_MM_ep.setTitleX("ep MM^2 (GeV^2)");
		H_dvcs_MM_eg = new H1F("H_dvcs_MM_eg","H_dvcs_MM_eg",100,-0.5,4);
		H_dvcs_MM_eg.setTitle("DVCS e#gamma MM^2");
		H_dvcs_MM_eg.setTitleX("e#gamma MM^2 (GeV^2)");
		H_dvcs_MM_epg = new H1F("H_dvcs_MM_epg","H_dvcs_MM_epg",100,-0.25,0.25);
		H_dvcs_MM_epg.setTitle("DVCS ep#gamma MM^2");
		H_dvcs_MM_epg.setTitleX("ep#gamma MM^2 (GeV^2)");
		H_dvcs_Phi = new H1F("H_dvcs_Phi","H_dvcs_Phi",100,-180,180);
		H_dvcs_Phi.setTitle("DVCS phi");
		H_dvcs_Phi.setTitleX("#phi (^o)");

		H_pi0_Q2xB = new H2F("H_pi0_Q2xB","H_pi0_Q2xB",100,0,1,100,0,12);
		H_pi0_Q2xB.setTitle("DV#pi^0 Q2 vs xB");
		H_pi0_Q2xB.setTitleX("xB");
		H_pi0_Q2xB.setTitleY("Q2");
		H_pi0_tphi = new H2F("H_pi0_tphi","H_pi0_tphi",100,-180,180,100,0,5);
		H_pi0_tphi.setTitle("DV#pi^0 -t vs #phi");
		H_pi0_tphi.setTitleX("#phi (^o)");
		H_pi0_tphi.setTitleY("-t (GeV)");
		H_gg_open_E = new H2F("H_gg_open_E","H_gg_open_E",100,0,14,100,0,20);
		H_gg_open_E.setTitle("DV#pi^0 #gamma#gamma open angle vs E");
		H_gg_open_E.setTitleX("E #pi^0 (GeV)");
		H_gg_open_E.setTitleY("#gamma#gamma open angle (^o)");
		H_pi0_pT = new H1F("H_pi0_pT","H_pi0_pT",100,0,2.5);
		H_pi0_pT.setTitle("DV#pi^0 missing pT");
		H_pi0_pT.setTitleX("pT (GeV)");
		H_pi0_cone = new H1F("H_pi0_cone","H_pi0_cone",100,0,25);
		H_pi0_cone.setTitle("DV#pi^0 cone angle");
		H_pi0_cone.setTitleX("cone ang (^o)");
		H_pi0_copl = new H1F("H_pi0_copl","H_pi0_copl",100,0,60);
		H_pi0_copl.setTitle("DV#pi^0 copl angle");
		H_pi0_copl.setTitleX("copl ang (^o)");
		H_pi0_ME = new H1F("H_pi0_ME","H_pi0_ME",100,-5,5);
		H_pi0_ME.setTitle("DV#pi^0 miss E");
		H_pi0_ME.setTitleX("miss E (GeV)");
		H_pi0_MM_ep = new H1F("H_pi0_MM_ep","H_pi0_MM_ep",100,-3.5,2.5);
		H_pi0_MM_ep.setTitle("DV#pi^0 ep MM^2");
		H_pi0_MM_ep.setTitleX("ep MM^2 (GeV^2)");
		H_pi0_MM_epi0 = new H1F("H_pi0_MM_epi0","H_pi0_MM_epi0",100,-2.5,7.5);
		H_pi0_MM_epi0.setTitle("DV#pi^0 e#pi^0 MM^2");
		H_pi0_MM_epi0.setTitleX("e#pi^0 MM^2 (GeV^2)");
		H_pi0_MM_eppi0 = new H1F("H_pi0_MM_eppi0","H_pi0_MM_eppi0",100,-2,1);
		H_pi0_MM_eppi0.setTitle("DV#pi^0 ep#pi^0 MM");
		H_pi0_MM_eppi0.setTitleX("ep#pi^0 MM^2 (GeV^2)");
		H_pi0_Phi = new H1F("H_pi0_Phi","H_pi0_Phi",100,-180,180);
		H_pi0_Phi.setTitle("DV#pi^0 #phi");
		H_pi0_Phi.setTitleX("#phi (^0)");
		H_IM_pi0 = new H1F("H_IM_pi0","H_IM_pi0",100,0,0.25);
		H_IM_pi0.setTitle("DV#pi^0 #gamma#gamma IM");
		H_IM_pi0.setTitleX("#gamma#gamma IM (GeV)");

		H_elast_e_th_mom = new H2F[7];
		H_elast_e_th_phi = new H2F[7];
		H_elast_e_phi_mom = new H2F[7];
		H_elast_W_th = new H2F[7];
		H_elast_W_phi = new H2F[7];
		H_elast_p_th_mom = new H2F[7];
		H_elast_p_th_phi = new H2F[7];
		H_elast_p_phi_mom = new H2F[7];
		H_elast_th_th = new H2F[7];
		H_elast_ph_ph = new H2F[7];
		H_elast_EB_bb = new H2F[7];
		H_elast_vz_vz = new H2F[7];
		for(int s=0;s<7;s++){
			H_elast_e_th_mom[s] = new H2F(String.format("H_elast_e_th_mom_S%d",s),String.format("S%d e #theta vs mom",s),100,0,EB,100,0,40);
			H_elast_e_th_phi[s] = new H2F(String.format("H_elast_e_th_phi_S%d",s),String.format("S%d e #theta vs #phi",s),100,-60,60,100,0,40);
			H_elast_e_phi_mom[s] = new H2F(String.format("H_elast_e_phi_mom_S%d",s),String.format("S%d e #phi vs mom",s),100,0,EB,100,-60,60);
			H_elast_W_th[s] = new H2F(String.format("H_elast_W_th_S%d",s),String.format("S%d W vs #theta",s),100,0,40,100,0,4);
			H_elast_W_phi[s] = new H2F(String.format("H_elast_W_phi_S%d",s),String.format("S%d W vs #phi",s),100,-60,60,100,0,4);
			H_elast_p_th_mom[s] = new H2F(String.format("H_elast_p_th_mom_S%d",s),String.format("S%d p #theta vs mom",s),100,0,4,100,0,70);
			H_elast_p_th_phi[s] = new H2F(String.format("H_elast_p_th_phi_S%d",s),String.format("S%d p #theta vs mom",s),100,-90,30,100,0,70);
			H_elast_p_phi_mom[s] = new H2F(String.format("H_elast_p_phi_mom_S%d",s),String.format("S%d p #phi vs mom",s),100,0,4,100,-90,30);
			H_elast_th_th[s] = new H2F(String.format("H_elast_th_th_S%d",s),String.format("S%d #theta vs #theta",s),100,0,70,100,0,40);
			H_elast_ph_ph[s] = new H2F(String.format("H_elast_ph_ph_S%d",s),String.format("S%d #phi vs #phi",s),100,-180,180,100,-180,180);
			H_elast_EB_bb[s] = new H2F(String.format("H_elast_EB_bb_S%d",s),String.format("S%d Ebeam vs #Delta#phi",s),100,-180,180,100,0,50);
			H_elast_vz_vz[s] = new H2F(String.format("H_elast_vz_vz_S%d",s),String.format("S%d vz p vs vz e",s),100,-25,25,100,-25,25);
		}
		H_dvcs_phi = new TH1F("H_dvcs_phi","H_dvcs_phi",NPHI,-180,180);
		H_dvcs_phi.setTitle("DVCS #phi counts");
		H_dvcs_phi.setTitleX("#phi (^o)");
		H_dvcs_phi_plus = new TH1F("H_dvcs_phi_plus","H_dvcs_phi_plus",NPHI,-180,180);
		H_dvcs_phi_plus.setTitle("DVCS #phi counts h+");
		H_dvcs_phi_plus.setTitleX("#phi (^o)");
		H_dvcs_phi_minus = new TH1F("H_dvcs_phi_minus","H_dvcs_phi_minus",NPHI,-180,180);
		H_dvcs_phi_minus.setTitle("DVCS #phi counts h-");
		H_dvcs_phi_minus.setTitleX("#phi (^o)");
		H_pi0_phi = new TH1F("H_pi0_phi","H_pi0_phi",NPHI,-180,180);
		H_pi0_phi.setTitle("#pi^0 #phi counts");
		H_pi0_phi.setTitleX("#phi (^o)");
		H_pi0_phi_plus = new TH1F("H_pi0_phi_plus","H_pi0_phi_plus",NPHI,-180,180);
		H_pi0_phi_plus.setTitle("#pi^0 #phi counts h+");
		H_pi0_phi_plus.setTitleX("#phi (^o)");
		H_pi0_phi_minus = new TH1F("H_pi0_phi_minus","H_pi0_phi_minus",NPHI,-180,180);
		H_pi0_phi_minus.setTitle("#pi^0 #phi counts h-");
		H_pi0_phi_minus.setTitleX("#phi (^o)");


		H_elec_mom = new TH1F("H_elec_mom", "H_elec_mom",100, 0, EB);
		H_elec_mom.setTitle("Electron momentum");
		H_elec_mom.setTitleX("p(#e) (GeV/c)");
		H_prot_mom = new TH1F("H_prot_mom", "H_prot_mom",100, 0, EB);
		H_prot_mom.setTitle("Proton momentum");
		H_prot_mom.setTitleX("p(#p) (GeV/c)");
		H_phot_mom = new TH1F("H_phot_mom", "H_phot_mom",100, 0.5, EB);
		H_phot_mom.setTitle("Photon momentum");
		H_phot_mom.setTitleX("p(#gamma) (GeV/c)");

		H_Q2 = new TH1F("H_Q2","H_Q2",100,0,12);
		H_Q2.setTitleX("Q2 (GeV^2)");
		H_xB = new TH1F("H_xB","H_xB",20,0,2);
		H_xB.setTitleX("xB");
		H_t = new TH1F("H_t","H_t",100,0,5);
		H_t.setTitleX("-t (GeV^2)");
		H_W = new TH1F("H_W","H_W",100,0,EB);
		H_W.setTitleX("W (GeV)");

		H_Q2xB = new TH2F("H_Q2xB","H_Q2xB",100,0,1,100,0,12);
		H_Q2xB.setTitle("Q2 vs xB");
		H_Q2xB.setTitleX("xB");
		H_Q2xB.setTitleY("Q2 (GeV^2)");
		H_tphi = new TH2F("H_tphi","H_tphi",100,-180,180,100,0,5);
		H_tphi.setTitle("-t vs #phi");
		H_tphi.setTitleX("#phi (^o)");
		H_tphi.setTitleY("-t (GeV^2)");

		H_MM_ep = new TH1F("H_MM_ep","H_MM_ep",100,-1,1.2);
		H_MM_ep.setTitle("ep MM^2");
		H_MM_ep.setTitleX("ep MM^2 (GeV^2)");
		H_MM_eg = new TH1F("H_MM_eg","H_MM_eg",100,-0.5,4);
		H_MM_eg.setTitle("e#gamma MM^2");
		H_MM_eg.setTitleX("e#gamma MM^2 (GeV^2)");
		H_MM_epg = new TH1F("H_MM_epg","H_MM_epg",100,-0.25,0.25);
		H_MM_epg.setTitle("ep#gamma MM^2");
		H_MM_epg.setTitleX("ep#gamma MM^2 (GeV^2)");

		H_elec_all_theta_phi = new TH2F("H_elec_all_theta_phi","H_elec_all_theta_phi",100,-180,180,100,0,40);
		H_elec_all_theta_phi.setTitle("elec #theta vs #phi");
		H_elec_all_theta_phi.setTitleX("#phi (^o)");
		H_elec_all_theta_phi.setTitleY("#theta (^o)");

		H_dvcs_elec_mom = new TH1F("H_dvcs_elec_mom", "H_dvcs_elec_mom",100, 0, EB);
		H_dvcs_elec_mom.setTitle("DVCS Electron momentum");
		H_dvcs_elec_mom.setTitleX("p(#e) (GeV/c)");
		H_dvcs_prot_mom = new TH1F("H_dvcs_prot_mom", "H_dvcs_prot_mom",100, 0, EB);
		H_dvcs_prot_mom.setTitle("DVCS Proton momentum");
		H_dvcs_prot_mom.setTitleX("p(#p) (GeV/c)");
		H_dvcs_phot_mom = new TH1F("H_dvcs_phot_mom", "H_dvcs_phot_mom",100, 0.5, EB);
		H_dvcs_phot_mom.setTitle("DVCS Photon momentum");
		H_dvcs_phot_mom.setTitleX("p(#gamma) (GeV/c)");

		H_dvcs_Q2 = new TH1F("H_dvcs_Q2","H_dvcs_Q2",100,0,12);
		H_dvcs_Q2.setTitleX("DVCS Q2 (GeV^2)");
		H_dvcs_xB = new TH1F("H_dvcs_xB","H_dvcs_xB",100,0,1);
		H_dvcs_xB.setTitleX("DVCS xB");
		H_dvcs_t = new TH1F("H_dvcs_t","H_dvcs_t",100,0,5);
		H_dvcs_t.setTitleX("DVCS -t (GeV^2)");
		H_dvcs_W = new TH1F("H_dvcs_W","H_dvcs_W",100,0,EB);
		H_dvcs_W.setTitleX("DVCS W (GeV)");

		H_dvcs_elec_theta_phi = new TH2F("H_dvcs_elec_theta_phi","H_dvcs_elec_theta_phi",100,-180,180,100,0,40);
		H_dvcs_elec_theta_phi.setTitle("DVCS elec #theta vs #phi");
		H_dvcs_elec_theta_phi.setTitleX("#phi (^o)");
		H_dvcs_elec_theta_phi.setTitleY("#theta (^o)");


}

void plot(){
		// -------------------------------------------------------------------------------------------
	// Drawing histograms
	TCanvas * c0 = new TCanvas();
	c0 -> Divide(2,1);
	c0 -> cd(1);
	h1_e_vz -> Draw(      );
	h1_p_vz -> Draw("same");
	c0 -> cd(2);
	h1_dlt_vz_ep   -> Draw(      );
	c0 -> Modified();
	c0 -> Update();

	TCanvas * c1 = new TCanvas();
	c1 -> Divide(2, 2);
	c1 -> cd(1);	h1_pmx  -> Draw();
	c1 -> cd(2);	h1_pmy  -> Draw();
	c1 -> cd(3);	h1_pmz  -> Draw();
	c1 -> cd(4);	h1_pmiss-> Draw();
	c1 -> Modified();
	c1 -> Update(); 

	TCanvas * c2 = new TCanvas();
	c2 -> Divide(2, 2);
	c2 -> cd(1);	h1_e_th     -> Draw(      );
	c2 -> cd(2);	h2_e_th_phi -> Draw("COLZ");
	c2 -> cd(4);	h1_e_phi    -> Draw(      );
	c2 -> Modified();
	c2 -> Update();

	TCanvas * c3 = new TCanvas();
	c3 -> Divide(2, 1);
	c3 -> cd(1);	h1_W  -> Draw();
	c3 -> cd(2);	h1_xB -> Draw();
	c3 -> Modified();
	c3 -> Update();

	TCanvas * c5 = new TCanvas();
	c5 -> Divide(2,1);
	c5 -> cd(1);	gPad -> SetLogz();	h2_e_Ep_p_0 -> Draw("COLZ");
	c5 -> cd(2);	gPad -> SetLogz();	h2_e_Ep_p_1 -> Draw("COLZ");
	c5 -> Modified();
	c5 -> Update();

	TCanvas * c6 = new TCanvas();
	c6 -> Divide(2, 2);
	c6 -> cd(1);	h1_e_px -> Draw();
	c6 -> cd(2);	h1_e_py -> Draw();
	c6 -> cd(3);	h1_e_pz -> Draw();
	c6 -> cd(4);	h1_e_p  -> Draw();
	c6 -> Modified();
	c6 -> Update();

	TCanvas * c7 = new TCanvas();
	c7 -> Divide(2, 2);
	c7 -> cd(1);	h1_p_px -> Draw();
	c7 -> cd(2);	h1_p_py -> Draw();
	c7 -> cd(3);	h1_p_pz -> Draw();
	c7 -> cd(4);	h1_p_p  -> Draw();
	c7 -> Modified();
	c7 -> Update();

	TCanvas * c8 = new TCanvas();
	c8->Divide(2, 1);
	c8->cd(1);	gPad -> SetLogz();	h2_beta_p_pos -> Draw("COLZ");
	c8->cd(2);	gPad -> SetLogz();	h2_beta_p_p   -> Draw("COLZ");
	c8 -> Modified();
	c8 -> Update(); 

	TCanvas * c9 = new TCanvas();
	c9 -> Divide(2, 2);
	c9 -> cd(1);	h1_p_th    -> Draw(      );
	c9 -> cd(2);	h2_p_th_phi-> Draw("COLZ");
	c9 -> cd(4);	h1_p_phi   -> Draw(      );
	c9 -> Modified();
	c9 -> Update();

	TCanvas * c10 = new TCanvas();
	c10 -> Divide(2, 2);
	c10 -> cd(1);	h1_e_lu -> Draw();
	c10 -> cd(2);	h1_e_lv -> Draw();
	c10 -> cd(3);	h1_e_lw -> Draw();
	c10 -> Modified();
	c10 -> Update();

	TCanvas * c11 = new TCanvas();
	h1_p_num -> Draw();
	c11 -> Modified();
	c11 -> Update();

	TCanvas * c12 = new TCanvas();
	h1_Mmiss -> Draw();
	c12 -> Modified();
	c12 -> Update();

	TCanvas * c13 = new TCanvas();
	c13 -> Divide(2,1);
	c13 -> cd(1);	h2_e_vz_phi-> Draw("COLZ");
	c13 -> cd(2);	h2_p_vz_phi-> Draw("COLZ");
	c13 -> Modified();
	c13 -> Update();

	TCanvas * c14 = new TCanvas();
	c14 -> Divide(2,1);
	c14 -> cd(1);	h1_e_tof  -> Draw(      );
	c14 -> cd(2);	h2_e_tof_p-> Draw("COLZ");
	c14 -> Modified();
	c14 -> Update();

	TCanvas * c15 = new TCanvas();
	c15 -> Divide(2, 1);
	c15 -> cd(1);	gPad -> SetLogz();	h2_p_dtT_p_0 -> Draw("COLZ");
	c15 -> cd(2);	gPad -> SetLogz();	h2_p_dtT_p_1 -> Draw("COLZ");
	c15 -> Modified();
	c15 -> Update();

	TCanvas * c16 = new TCanvas();
	h2_p_tof_det -> Draw("COLZ");
	c16 -> Modified();
	c16 -> Update();

	TCanvas * c17 = new TCanvas();
	h2_p_dtT_det -> Draw("COLZ");
	c17 -> Modified();
	c17 -> Update();

	TCanvas * c18 = new TCanvas();
	h1_Em -> Draw();
	c18 -> Modified();
	c18 -> Update();

	TCanvas * c19 = new TCanvas();
	h2_Em_Pm -> Draw("COLZ");
	c19 -> Modified();
	c19 -> Update();

	TCanvas * c20 = new TCanvas();
	h2_pe_pp -> Draw("COLZ");
	c20 -> Modified();
	c20 -> Update();

	TCanvas * c21 = new TCanvas();
	h1_pm_th -> Draw();
	c21 -> Modified();
	c21 -> Update();

	TCanvas * c22 = new TCanvas();
	h1_pm_ph -> Draw();
	c22 -> Modified();
	c22 -> Update();

	TCanvas * c23 = new TCanvas();
	h2_pe_the -> Draw("COLZ");
	c23 -> Modified();
	c23 -> Update();

	TCanvas * c24 = new TCanvas();
	c24 -> Divide(2,1);
	c24 -> cd(1);	h1_p_th_meas_calc -> Draw();
	c24 -> cd(2);	h1_p_p_meas_calc  -> Draw();
	c24 -> Modified();
	c24 -> Update();

	// -------------------------------------------------------------------------------------------
	// Saving plots to system

	c0  -> Print("results_eP.pdf(");
	c1  -> Print("results_eP.pdf" );
	c2  -> Print("results_eP.pdf" );
	c3  -> Print("results_eP.pdf" );
	c5  -> Print("results_eP.pdf" );
	c6  -> Print("results_eP.pdf" );
	c7  -> Print("results_eP.pdf" );
	c8  -> Print("results_eP.pdf" );
	c9  -> Print("results_eP.pdf" );
	c10 -> Print("results_eP.pdf" );
	c11 -> Print("results_eP.pdf" );
	c12 -> Print("results_eP.pdf" );
	c13 -> Print("results_eP.pdf" );
	c14 -> Print("results_eP.pdf" );
	c15 -> Print("results_eP.pdf" );
	c16 -> Print("results_eP.pdf" );
	c17 -> Print("results_eP.pdf" );
	c18 -> Print("results_eP.pdf" );
	c19 -> Print("results_eP.pdf" );
	c20 -> Print("results_eP.pdf" );
	c21 -> Print("results_eP.pdf" );
	c22 -> Print("results_eP.pdf" );
	c23 -> Print("results_eP.pdf" );
	c24 -> Print("results_eP.pdf)");

}