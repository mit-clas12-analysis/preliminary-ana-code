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


bool isMC, isDVCS, isPI0, usePCAL, useFTCAL;
// Random myMC;
double EB, HWP;
int HELI, NPHI, haz_elec, haz_prot, haz_g1, haz_g2, runNum, Ndvcs, Npi0;
int e_sect, p_sect;
double elec_mom, elec_theta, elec_phi, elec_vz, elec_phi_sect;
double prot_mom, prot_theta, prot_phi, prot_vz, prot_beta, prot_phi_sect;
double g1_E, g1_theta, g1_phi, g1_beta;
double g2_E, g2_theta, g2_phi, g2_beta;
double elast_W, elast_EB, elast_bb, all_W;

TLorentzVector VB, VT, VE, VGS, VPROT, Vmand, VG1, VG2, VPI0, VMISS, VmissP, VmissG;
TVector3 Vlept, Vhadr, Vhad2;

TH2F **H_elec_theta_mom, **H_elec_theta_phi, **H_elec_phi_mom;
TH2F **H_prot_theta_mom, **H_prot_theta_phi, **H_prot_phi_mom;
TH2F **H_g_theta_E, **H_g_theta_phi, **H_g_phi_E;
TH2F **H_g1_theta_E, **H_g1_theta_phi, *H_g1_phi_E;
TH2F *H_g2_theta_E, *H_g2_theta_phi, *H_g2_phi_E;
TH2F *H_prot_beta_p, *H_prot_e_vz, *H_phot_b_p;

TH2F *H_prot_miss_phi, *H_phot_elec_angle, *H_phot_virt_angle;
TH2F *H_prot_pi0_miss_phi, *H_pi0_elec_angle, *H_pi0_virt_angle;

TH2F *H_dvcs_Q2xB, *H_dvcs_tphi;
TH1F *H_dvcs_pT, *H_dvcs_cone, *H_dvcs_copl, *H_dvcs_ME, *H_dvcs_MM_ep, *H_dvcs_MM_eg, *H_dvcs_MM_epg, *H_dvcs_Phi;
TH2F *H_pi0_Q2xB, *H_pi0_tphi, *H_gg_open_E;
TH1F *H_pi0_pT, *H_pi0_cone, *H_pi0_copl, *H_pi0_ME, *H_pi0_MM_ep, *H_pi0_MM_epi0, *H_pi0_MM_eppi0, *H_pi0_Phi, *H_IM_pi0;
	GraphErrors g_dvcs_bsa;

TH2F *H_elast_e_th_mom[7], *H_elast_e_th_phi[7], *H_elast_e_phi_mom[7];
TH2F *H_elast_p_th_mom[7], *H_elast_p_th_phi[7], *H_elast_p_phi_mom[7];
TH2F *H_elast_vz_vz[7], *H_elast_W_th[7], *H_elast_W_phi[7], *H_elast_th_th[7], *H_elast_ph_ph[7], *H_elast_EB_bb[7];

TH1F *H_dvcs_phi, *H_dvcs_phi_plus, *H_dvcs_phi_minus;
TH1F *H_pi0_phi, *H_pi0_phi_plus, *H_pi0_phi_minus;

	// electron mom, proton mom, photon mom
	// Q2 xb W, t
	// Q2 vs W, t vs phi
	// angle variables theta phi
	// missing mass ep epgamma
	// entire, dvcs cut

TH1F *H_elec_mom, *H_prot_mom, *H_phot_mom;
TH1F *H_Q2, *H_xB, *H_W, *H_t;
TH2F *H_Q2xB, *H_tphi;
TH1F *H_MM_ep, *H_MM_eg, *H_MM_epg;
TH2F *H_elec_all_theta_phi;

TH1F *H_dvcs_elec_mom, *H_dvcs_prot_mom, *H_dvcs_phot_mom;
TH1F *H_dvcs_Q2, *H_dvcs_xB, *H_dvcs_W, *H_dvcs_t;
TH2F *H_dvcs_elec_theta_phi;
