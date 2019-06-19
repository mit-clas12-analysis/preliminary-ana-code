#include <cstdlib>
#include <iostream>
#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "reader.h"
#include "dictionary.h"
#include "particle.h"
#include "calorimeter.h"
#include "scintillator.h"
#include "band.h"
#include "event.h"

using namespace std;

void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);
void InitiateHistograms();
void plot();
TH1F * h1_e_vz;
TH1F * h1_e_tof;
TH1F * h1_e_px;
TH1F * h1_e_py;
TH1F * h1_e_pz;
TH1F * h1_e_p ;
TH1F * h1_e_th;
TH1F * h1_e_phi; 
TH1F * h1_e_lu;
TH1F * h1_e_lv;
TH1F * h1_e_lw;
TH1F * h1_p_vz;
TH1F * h1_p_num;
TH1F * h1_p_px;
TH1F * h1_p_py;
TH1F * h1_p_pz;
TH1F * h1_p_p ;
TH1F * h1_p_th;
TH1F * h1_p_phi ;
TH1F * h1_pmiss ;
TH1F * h1_pmx ;
TH1F * h1_pmy ;
TH1F * h1_pmz ;
TH1F * h1_pm_th ;
TH1F * h1_pm_ph ;
TH1F * h1_Mmiss ;
TH1F * h1_Em  ;
TH1F * h1_W   ;
TH1F * h1_xB  ;
TH1F * h1_dlt_vz_ep;
TH1F * h1_p_th_meas_calc;
TH1F * h1_p_p_meas_calc ;
// 2D histograms
TH2F * h2_e_Ep_p_0 ;
TH2F * h2_e_Ep_p_1 ;
TH2F * h2_e_th_phi ;
TH2F * h2_p_th_phi ;
TH2F * h2_beta_p_pos ;
TH2F * h2_beta_p_p ;
TH2F * h2_e_vz_phi ;
TH2F * h2_p_vz_phi ;
TH2F * h2_e_tof_p  ;
TH2F * h2_p_dtT_p_0;
TH2F * h2_p_dtT_p_1;
TH2F * h2_p_tof_det;
TH2F * h2_p_dtT_det;
TH2F * h2_Em_Pm    ;
TH2F * h2_pe_pp    ;
TH2F * h2_pe_the   ;