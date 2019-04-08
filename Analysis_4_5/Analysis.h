#ifndef zeeAna_h
#define zeeAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include "TH1F.h"
#include <iostream>
#include <vector>
#include "makeHists.h"
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TH3F.h"
#include <TRandom3.h>
#include <TMinuit.h>
#include <TApplication.h>
#include "TEnv.h"
#include <TComplex.h>
#include "TH2D.h"

using namespace std;
//class: the main class for functions;
class Analysis
{
 public:

  TTree      *fChain;
//*******variables for tree*******//
    //mc global
    

  Int_t           run;
  Int_t           lumi;
  ULong64_t       evt;
  Int_t           isData;
  Float_t         evt_scale1fb;
  Float_t         genps_weight;
  Float_t         xsec_br;
  Int_t           evt_passgoodrunlist;
  TString         *CMS4path;
  Int_t           CMS4index;
  Float_t         weight_fr_r1_f1;
  Float_t         weight_fr_r1_f2;
  Float_t         weight_fr_r1_f0p5;
  Float_t         weight_fr_r2_f1;
  Float_t         weight_fr_r2_f2;
  Float_t         weight_fr_r2_f0p5;
  Float_t         weight_fr_r0p5_f1;
  Float_t         weight_fr_r0p5_f2;
  Float_t         weight_fr_r0p5_f0p5;
  Float_t         weight_pdf_up;
  Float_t         weight_pdf_down;
  Float_t         weight_alphas_down;
  Float_t         weight_alphas_up;
  Int_t           HLT_DoubleMu;
  Int_t           HLT_DoubleEl;
  Int_t           HLT_MuEG;
  Int_t           pass_duplicate_ee_em_mm;
  Int_t           pass_duplicate_mm_em_ee;
  Float_t         gen_ht;
  Int_t           firstgoodvertex;
  Int_t           nvtx;
  Int_t           nTrueInt;
//  Int_t           lep_p4_;
//  Float_t         lep_p4_fCoordinates_fX[20];   //[lep_p4_]
//  Float_t         lep_p4_fCoordinates_fY[20];   //[lep_p4_]
//  Float_t         lep_p4_fCoordinates_fZ[20];   //[lep_p4_]
//  Float_t         lep_p4_fCoordinates_fT[20];   //[lep_p4_]
  vector<float>   *lep_pt=0;
  vector<float>   *lep_eta=0;
  vector<float>   *lep_phi=0;
  vector<float>   *lep_energy=0;
  vector<float>   *lep_mva=0;
  vector<float>   *lep_relIso03EA=0;
  vector<float>   *lep_relIso03EAwLep=0;
  vector<float>   *lep_ip3d=0;
  vector<float>   *lep_sip3d=0;
  vector<float>   *lep_dxy=0;
  vector<float>   *lep_dz=0;
  vector<int>     *lep_mc_id=0;
  vector<int>     *lep_motherIdv2=0;
  vector<int>     *lep_idx=0;
  vector<int>     *lep_id=0;
  vector<int>     *lep_isTightPOG=0;
  vector<int>     *lep_isMediumPOG=0;

//  Float_t         fCoordinates_fX;
//  Float_t         fCoordinates_fY;
//  Float_t         fCoordinates_fZ;
//  Float_t         fCoordinates_fT;
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_up_pt;
  Float_t         met_up_phi;
  Float_t         met_dn_pt;
  Float_t         met_dn_phi;
  Float_t         met_gen_pt;
  Float_t         met_gen_phi;
//  Int_t           jets_p4_;
//  Float_t         jets_p4_fCoordinates_fX[20];   //[jets_p4_]
//  Float_t         jets_p4_fCoordinates_fY[20];   //[jets_p4_]
//  Float_t         jets_p4_fCoordinates_fZ[20];   //[jets_p4_]
//  Float_t         jets_p4_fCoordinates_fT[20];   //[jets_p4_]
  vector<float>   *jets_pt=0;
  vector<float>   *jets_eta=0;
  vector<float>   *jets_phi=0;
  vector<float>   *jets_mass=0;
  Int_t           nj;
  Int_t           nb;
  Int_t           nbmed;
  Float_t         ht;
  Float_t         weight_btagsf;
  Float_t         weight_btagsf_heavy_DN;
  Float_t         weight_btagsf_heavy_UP;
  Float_t         weight_btagsf_light_DN;
  Float_t         weight_btagsf_light_UP;


  Float_t         lumi_weight = 41.3;

//*******Statistic variables******//
/*    double weighted_number;
*/
//*******control variables*******//
    bool isRead;
    bool isDamaged;

//*******DQ bad runs******//
    vector<int> DQBadRunList;

//******initial plots*******//

      TH1D *puw_2017_central;
      TH2D *SFP_muon_recoid;
      TH2D *SFP_elec_recoid;
      TH2D *SFP_elec_mva;
      TH2D *SFP_elec_mva_3l;
      TH2D *SFP_muon_isoip;
      TH2D *SFP_muon_isoip_3l;
      TH2D *SFP_elec_isoip;
      TH2D *SFP_elec_isoip_3l;
/*      TH2D *histmap_lead_mu_recoid_sf;
      TH2D *histmap_subl_mu_recoid_sf;
      TH2D *histmap_lead_el_recoid_sf;
      TH2D *histmap_subl_el_recoid_sf;
      TH2D *histmap_lead_el_mva_sf;
      TH2D *histmap_subl_el_mva_sf;
      TH2D *histmap_emu_mu_recoid_sf;
      TH2D *histmap_emu_el_recoid_sf;
      TH2D *histmap_emu_el_mva_sf;
      TH2D *histmap_lead_mu_recoid_3l_sf;
      TH2D *histmap_subl_mu_recoid_3l_sf;
      TH2D *histmap_lead_el_recoid_3l_sf;
      TH2D *histmap_subl_el_recoid_3l_sf;
      TH2D *histmap_lead_el_mva_3l_sf;
      TH2D *histmap_subl_el_mva_3l_sf;
      TH2D *histmap_tert_mu_recoid_3l_sf;
      TH2D *histmap_tert_el_recoid_3l_sf;
      TH2D *histmap_tert_el_mva_3l_sf;
      TH2D *histmap_lead_mu_isoip_sf;
      TH2D *histmap_subl_mu_isoip_sf;
      TH2D *histmap_lead_el_isoip_sf;
      TH2D *histmap_subl_el_isoip_sf;
      TH2D *histmap_emu_mu_isoip_sf;
      TH2D *histmap_emu_el_isoip_sf;
      TH2D *histmap_lead_mu_isoip_3l_sf;
      TH2D *histmap_subl_mu_isoip_3l_sf;
      TH2D *histmap_lead_el_isoip_3l_sf;
      TH2D *histmap_subl_el_isoip_3l_sf;
      TH2D *histmap_tert_mu_isoip_3l_sf;
      TH2D *histmap_tert_el_isoip_3l_sf;*/
//    TH1D *_Lumi_new_file;
//    TH2D *Zpt_2D_rwt;

//******effciency plots******//
//    TH2D *eff_muon_isomtc_pt_eta;
//    TH2D *eff_muon_trkmatch_pt_eta;
//    TH1D *eff_muon_iso;
//    TH1D *eff_muon_id;
//    TH1D *eff_muon_trk;
//    TH1D *eff_muon_trkmatch_eta;
//    TH2D *eff_muon_isomtc_eta_trkz;
//    TH2D *eff_muon_trkmatch_eta_trkz;
//    TH2D *eff_muon_isomtc_eta_phi;
//    TH2D *trigger_eff_etaphibin;

//*******random number*****//
//    TRandom3 *MyR_GS;
//*******functions********//
  Analysis(const char* ifileName, const char* TypeName);
  virtual ~Analysis();
  virtual void  Initial(const char* RootName, int RootNumber);
  virtual void  Loop(const char* TypeName);
  virtual void  End(int RootNumber);
  virtual void  Finish(int RootNumber);
  virtual void  Output();

//*******plots*******//
  makeHists* myhists;

//*******test and output files*******//

//  ofstream outfile;
//  outfile.open("test.txt",ios::out);

};
#endif

#ifdef zeeAna_cxx
Analysis::Analysis(const char* ifileName, const char* TypeName)
{
//initializing random
// MyR_GS = new TRandom3(54321);

//initializing histograms
  myhists = new makeHists();
  TString histoname = (TString)TypeName + "_results.root";
  myhists->bookHists(histoname);

//initializing statistics
/*  weighted_number = 0;
  Sta_InputNumber = 0;
  Sta_PreNumber_PassDQ = 0;
  Sta_PreNumber_FailDQ = 0;
  Sta_PreNumber_total = 0;
  Sta_passNumber_PassDQ = 0;
  Sta_passNumber_FailDQ = 0;
  Sta_passNumber_total = 0;
  Sta_SSevt = 0;
  Sta_OSevt = 0;
  Sta_SSevt_TypeA = 0;
  Sta_OSevt_TypeA = 0;
  Sta_SSevt_TypeB = 0;
  Sta_OSevt_TypeB = 0;
  Sta_SSevt_TypeC = 0;
  Sta_OSevt_TypeC = 0;
*/
/*  for(int icmp=0;icmp<4;icmp++)
   {
    Sta_SSevt_Eta[icmp]=0;
    Sta_OSevt_Eta[icmp]=0;
   }

  for(int iSS=0;iSS<8;iSS++)
   {
    Sta_PtSMT[iSS]=0;
   }
*/
//initializing controls

//initializing global plots
  //initializing luminosity
/*   std::string lumi_new = "Inputfiles/Reweight/Lum_reweight_iia.C";                                      //for test
   std::string temp2 = string(".x ") + lumi_new;
   gROOT->ProcessLine(temp2.c_str());
   _Lumi_new_file = (TH1D *)gROOT->FindObject("Scalefactor");
   _Lumi_new_file->SetName("lumi reweight");
*/
  
   TFile *puw_central = new TFile("./scalefactors/puw_2017.root","READ");
   puw_2017_central = (TH1D*)(puw_central->Get("puw_central"));

   TFile *FMReco     = new TFile("./scalefactors/RunBCDEF_SF_ID.root","READ");
   SFP_muon_recoid   = (TH2D*)(FMReco->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   TFile *FEReco     = new TFile("./scalefactors/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","READ");
   SFP_elec_recoid   = (TH2D*)(FEReco->Get("EGamma_SF2D"));
   TFile *FEMva      = new TFile("./scalefactors/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp80noiso.root","READ");
   SFP_elec_mva      = (TH2D*)(FEMva->Get("EGamma_SF2D"));
   TFile *FEMva3l    = new TFile("./scalefactors/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp90noiso.root","READ");
   SFP_elec_mva_3l   = (TH2D*)(FEMva3l->Get("EGamma_SF2D"));
   TFile *FMIso      = new TFile("./scalefactors/isoipsf/MuonID_2017www/muon/MuMediumPOG_MuTightVVV/sf.root","READ");
   SFP_muon_isoip    = (TH2D*)(FMIso->Get("h_sf_pt_vs_eta"));
   TFile *FMIso3l    = new TFile("./scalefactors/isoipsf/MuonID_2017www/muon/MuMediumPOG_MuTightVVV3l/sf.root","READ");
   SFP_muon_isoip_3l = (TH2D*)(FMIso3l->Get("h_sf_pt_vs_eta"));
   TFile *FEIso      = new TFile("./scalefactors/isoipsf/ElectronID_2017www/electron/EGammaMVA80POG2017_EGammaTightVVV/sf.root","READ");
   SFP_elec_isoip    = (TH2D*)(FEIso->Get("h_sf_pt_vs_eta"));
   TFile *FEIso3l    = new TFile("./scalefactors/isoipsf/ElectronID_2017www/electron/EGammaMVA90POG2017_EGammaTightVVV3l/sf.root","READ");
   SFP_elec_isoip_3l = (TH2D*)(FEIso3l->Get("h_sf_pt_vs_eta"));

/*
   TFile *FRunBCDEF_SF_ID = new TFile("./scalefactors/RunBCDEF_SF_ID.root","READ");
   histmap_lead_mu_recoid_sf    = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   histmap_subl_mu_recoid_sf    = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   histmap_emu_mu_recoid_sf     = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   histmap_lead_mu_recoid_3l_sf = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   histmap_subl_mu_recoid_3l_sf = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
   histmap_tert_mu_recoid_3l_sf = (TH2D*)(FRunBCDEF_SF_ID->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));

   TFile *FegammaEffRECO = new TFile("./scalefactors/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root","READ");
   histmap_lead_el_recoid_sf    = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));
   histmap_subl_el_recoid_sf    = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));
   histmap_emu_el_recoid_sf     = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));
   histmap_lead_el_recoid_3l_sf = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));
   histmap_subl_el_recoid_3l_sf = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));
   histmap_tert_el_recoid_3l_sf = (TH2D*)(FegammaEffRECO->Get("EGamma_SF2D"));

   TFile *FgammaEff80 = new TFile("./scalefactors/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp80noiso.root","READ");
   histmap_lead_el_mva_sf = (TH2D*)(FgammaEff80->Get("EGamma_SF2D"));
   histmap_subl_el_mva_sf = (TH2D*)(FgammaEff80->Get("EGamma_SF2D"));
   histmap_emu_el_mva_sf  = (TH2D*)(FgammaEff80->Get("EGamma_SF2D"));

   TFile *FgammaEff90 = new TFile("./scalefactors/gammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp90noiso.root","READ");
   histmap_lead_el_mva_3l_sf = (TH2D*)(FgammaEff90->Get("EGamma_SF2D"));
   histmap_subl_el_mva_3l_sf = (TH2D*)(FgammaEff90->Get("EGamma_SF2D"));
   histmap_tert_el_mva_3l_sf = (TH2D*)(FgammaEff90->Get("EGamma_SF2D"));

   TFile *FMuMedium = new TFile("./scalefactors/isoipsf/MuonID_2017www/muon/MuMediumPOG_MuTightVVV/sf.root","READ");
   histmap_lead_mu_isoip_sf    = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));
   histmap_subl_mu_isoip_sf    = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));
   histmap_emu_mu_isoip_sf     = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));
   histmap_lead_mu_isoip_3l_sf = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));
   histmap_subl_mu_isoip_3l_sf = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));
   histmap_tert_mu_isoip_3l_sf = (TH2D*)(FMuMedium->Get("h_sf_pt_vs_eta"));

   TFile *FEGammaVVV = new TFile("./scalefactors/isoipsf/ElectronID_2017www/electron/EGammaMVA80POG2017_EGammaTightVVV/sf.root","READ");
   histmap_lead_el_isoip_sf = (TH2D*)(FEGammaVVV->Get("h_sf_pt_vs_eta"));
   histmap_subl_el_isoip_sf = (TH2D*)(FEGammaVVV->Get("h_sf_pt_vs_eta"));
   histmap_emu_el_isoip_sf  = (TH2D*)(FEGammaVVV->Get("h_sf_pt_vs_eta"));

   TFile *FEGammaVVV3l = new TFile("./scalefactors/isoipsf/ElectronID_2017www/electron/EGammaMVA90POG2017_EGammaTightVVV3l/sf.root","READ");
   histmap_lead_el_isoip_3l_sf = (TH2D*)(FEGammaVVV3l->Get("h_sf_pt_vs_eta"));
   histmap_subl_el_isoip_3l_sf = (TH2D*)(FEGammaVVV3l->Get("h_sf_pt_vs_eta"));
   histmap_tert_el_isoip_3l_sf = (TH2D*)(FEGammaVVV3l->Get("h_sf_pt_vs_eta"));

*/

//   TFile *FSFLMRe = new TFile("./scalefactors/")

  //2D zpt reweighting
//   string temp2 = ".x Inputfiles/Reweight/Resbos_Pythia_ZPt_Y_Ratio.C";
//   gROOT->ProcessLine(temp2.c_str());
//   Zpt_2D_rwt = (TH2D *)gROOT->FindObject("_h_zypt");

//   TFile *lumi_file = new TFile("luminosity_iia.root","READ");
//   _Lumi_new_file = (TH1D*)(lumi_file->Get("lumi_scale"));

//   TFile *eff_scale = new TFile("ScaleFactor.root","READ");
//   eff_muon_id = (TH1D*)(eff_scale->Get("IdScale"));
//   eff_muon_iso = (TH1D*)(eff_scale->Get("IsoScale"));
//   eff_muon_trk = (TH1D*)(eff_scale->Get("TrkScale"));

/*
   string temp2 = ".x Inputfiles/Efficiency_iia/Id_Efficiency_Eta.C";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_id= (TH1D *)gROOT->FindObject("Id_Eff_Depend_Eta");
   eff_muon_id->SetName("Id_Eff");


   temp2 = ".x Inputfiles/Efficiency_iia/Iso_Efficiency_Eta.C";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_iso  = (TH1D *)gROOT->FindObject("Iso_Eff_Depend_Eta");
   eff_muon_iso->SetName("Iso_Eff");


   temp2 = ".x Inputfiles/Efficiency_iia/Trk_Efficiency_Eta.C";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_trk  = (TH1D *)gROOT->FindObject("Trk_Eff_Depend_Eta");
   eff_muon_trk->SetName("Trk_Eff");
*/
/*  //initializing efficiency
   temp2 = ".x Inputfiles/Efficiency/eff_muon_isomtc_pt_eta.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_isomtc_pt_eta = (TH2D *)gROOT->FindObject("eff_muon_isomtc_pt_eta");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_trkmatch_pt_eta.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_trkmatch_pt_eta = (TH2D *)gROOT->FindObject("eff_muon_trkmatch_pt_eta");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_isomtc_eta_trkz.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_isomtc_eta_trkz = (TH2D *)gROOT->FindObject("eff_muon_isomtc_eta_trkz");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_trkmatch_eta_trkz.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_trkmatch_eta_trkz = (TH2D *)gROOT->FindObject("eff_muon_trkmatch_eta_trkz");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_isomtc_eta_phi.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_isomtc_eta_phi = (TH2D *)gROOT->FindObject("eff_muon_isomtc_eta_phi");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_trkmatch_eta.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_trkmatch_eta = (TH1D *)gROOT->FindObject("scale_muon_trakmatch_eta");

   temp2 = ".x Inputfiles/Efficiency/eff_muon_isomtc_eta.cpp";
   gROOT->ProcessLine(temp2.c_str());
   eff_muon_isomtc_eta = (TH1D *)gROOT->FindObject("scale_muon_isomtc_eta");
*/                  
//   temp2 = ".x Inputfiles/Efficiency/trigger_eff_etaphibin.cpp";                                     
//   gROOT->ProcessLine(temp2.c_str());
//   trigger_eff_etaphibin = (TH2D *)gROOT->FindObject("eff_etaphibin");

//initializing DQ bad run list
/*  DQBadRunList.clear();
  ifstream DQinfile;
   DQinfile.open("Inputfiles/DQbadRuns/all-b2.runs");
  int getnumber = 0;
  while(DQinfile>>getnumber)
   {
    DQBadRunList.push_back(getnumber);
   }

  cout<<"DQ bad run list:   "<<DQBadRunList.size()<<endl;
*/
}

void Analysis::End(int RootNumber)
{
//output: free the current file
 cout<<"**Runing: Free Rootfile "<<RootNumber<<endl; 

 if(!fChain) 
  {
//output: reading file failed
   cout<<"XXXX**Runing: Big Error! No file load!"<<endl;
   return;
  }
 delete fChain->GetCurrentFile();

}


Analysis::~Analysis()
{
}

void Analysis::Finish(int RootNumber)
{
 //for AFB plots
/* for(int ibin=0;ibin<8;ibin++)
  {
   double bin_forward = myhists->EvtA_select_Zboson_mass_forward->GetBinContent(ibin+1);
   double bin_backward = myhists->EvtA_select_Zboson_mass_backward->GetBinContent(ibin+1);
   double err_forward = myhists->EvtA_select_Zboson_mass_forward->GetBinError(ibin+1);
   double err_backward = myhists->EvtA_select_Zboson_mass_backward->GetBinError(ibin+1);
   double set_bin = (bin_forward-bin_backward)/(bin_forward + bin_backward);
   double set_err = sqrt( (4.0*bin_backward*bin_backward / pow(bin_forward+bin_backward ,4))*err_forward*err_forward
                      + (4.0*bin_forward*bin_forward / pow(bin_forward+bin_backward ,4))*err_backward*err_backward);

   myhists->EvtA_select_AFB_mass->SetBinContent(ibin+1, set_bin);
   myhists->EvtA_select_AFB_mass->SetBinError(ibin+1, set_err);

   bin_forward = myhists->EvtB_select_Zboson_mass_forward->GetBinContent(ibin+1);
   bin_backward = myhists->EvtB_select_Zboson_mass_backward->GetBinContent(ibin+1);
   err_forward = myhists->EvtB_select_Zboson_mass_forward->GetBinError(ibin+1);
   err_backward = myhists->EvtB_select_Zboson_mass_backward->GetBinError(ibin+1);
   set_bin = (bin_forward-bin_backward)/(bin_forward + bin_backward);
   set_err = sqrt( (4.0*bin_backward*bin_backward / pow(bin_forward+bin_backward ,4))*err_forward*err_forward
                      + (4.0*bin_forward*bin_forward / pow(bin_forward+bin_backward ,4))*err_backward*err_backward);

   myhists->EvtB_select_AFB_mass->SetBinContent(ibin+1, set_bin);
   myhists->EvtB_select_AFB_mass->SetBinError(ibin+1, set_err);

   bin_forward = myhists->EvtC_select_Zboson_mass_forward->GetBinContent(ibin+1);
   bin_backward = myhists->EvtC_select_Zboson_mass_backward->GetBinContent(ibin+1);
   err_forward = myhists->EvtC_select_Zboson_mass_forward->GetBinError(ibin+1);
   err_backward = myhists->EvtC_select_Zboson_mass_backward->GetBinError(ibin+1);
   set_bin = (bin_forward-bin_backward)/(bin_forward + bin_backward);
   set_err = sqrt( (4.0*bin_backward*bin_backward / pow(bin_forward+bin_backward ,4))*err_forward*err_forward
                      + (4.0*bin_forward*bin_forward / pow(bin_forward+bin_backward ,4))*err_backward*err_backward);

   myhists->EvtC_select_AFB_mass->SetBinContent(ibin+1, set_bin);
   myhists->EvtC_select_AFB_mass->SetBinError(ibin+1, set_err);

   bin_forward = myhists->EvtTotal_select_Zboson_mass_forward->GetBinContent(ibin+1);
   bin_backward = myhists->EvtTotal_select_Zboson_mass_backward->GetBinContent(ibin+1);
   err_forward = myhists->EvtTotal_select_Zboson_mass_forward->GetBinError(ibin+1);
   err_backward = myhists->EvtTotal_select_Zboson_mass_backward->GetBinError(ibin+1);
   set_bin = (bin_forward-bin_backward)/(bin_forward + bin_backward);
   set_err = sqrt( (4.0*bin_backward*bin_backward / pow(bin_forward+bin_backward ,4))*err_forward*err_forward
                      + (4.0*bin_forward*bin_forward / pow(bin_forward+bin_backward ,4))*err_backward*err_backward);

   myhists->EvtTotal_select_AFB_mass->SetBinContent(ibin+1, set_bin);
   myhists->EvtTotal_select_AFB_mass->SetBinError(ibin+1, set_err);
  }*/

//save histogram
 cout<<"**Runing: "<<RootNumber<<" rootfiles finished"<<endl;
 if(myhists)
  {
   myhists->saveHists();
   delete myhists;
  }
 cout<<"**Runing: histogram saved"<<endl;

}

void Analysis::Output()
{
//output: final results
/* double misCC,misEC;
 cout<<"**Runing: output the statistical information:"<<endl;
 cout<<"Total event input:            "<<Sta_InputNumber<<endl;
 cout<<"weighted events:              "<<weighted_number<<endl;
 cout<<"normalized factor:            "<<weighted_number/(double)Sta_InputNumber<<endl;
 cout<<"pre total number:             "<<Sta_PreNumber_total<<endl; 
 cout<<"pre pass DQ number:           "<<Sta_PreNumber_PassDQ<<endl; 
 cout<<"pre fail DQ number:           "<<Sta_PreNumber_FailDQ<<endl;
 cout<<"selected total number:        "<<Sta_passNumber_total<<endl;
 cout<<"selected pass DQ number:      "<<Sta_passNumber_PassDQ<<endl;
 cout<<"selected fail DQ number:      "<<Sta_passNumber_FailDQ<<endl;
 cout<<"same sign event:              "<<Sta_SSevt<<endl;
 cout<<"opposite sign event:          "<<Sta_OSevt<<endl;
 cout<<"pT SMT distribution:          "<<endl;
 cout<<"      pT<25, with SMT         "<<Sta_PtSMT[0]<<endl;
 cout<<"      pT<35, with SMT         "<<Sta_PtSMT[1]<<endl;
 cout<<"      pT<55, with SMT         "<<Sta_PtSMT[2]<<endl;
 cout<<"      pT>55, with SMT         "<<Sta_PtSMT[3]<<endl;
 cout<<"      pT<25, no SMT           "<<Sta_PtSMT[4]<<endl;
 cout<<"      pT<35, no SMT           "<<Sta_PtSMT[5]<<endl;
 cout<<"      pT<55, no SMT           "<<Sta_PtSMT[6]<<endl;
 cout<<"      pT>55, no SMT           "<<Sta_PtSMT[7]<<endl;

 cout<<"**************For All Kinds of Events**********"<<endl;
 cout<<"SS event type A:              "<<Sta_SSevt_TypeA<<endl;
 cout<<"OS event type A:              "<<Sta_OSevt_TypeA<<endl;
 cout<<"SS event type B:              "<<Sta_SSevt_TypeB<<endl;
 cout<<"OS event type B:              "<<Sta_OSevt_TypeB<<endl;
 cout<<"SS event type C:              "<<Sta_SSevt_TypeC<<endl;
 cout<<"OS event type C:              "<<Sta_OSevt_TypeC<<endl;

 cout<<"SS event |eta|<1.0:           "<<Sta_SSevt_Eta[0]<<endl;
 cout<<"OS event |eta|<1.0:           "<<Sta_OSevt_Eta[0]<<endl;
 cout<<"SS event 1.0<|eta|<1.6:       "<<Sta_SSevt_Eta[1]<<endl;
 cout<<"OS event 1.0<|eta|<1.6:       "<<Sta_OSevt_Eta[1]<<endl;
 cout<<"SS event |eta|>1.6:           "<<Sta_SSevt_Eta[2]<<endl;
 cout<<"OS event |eta|>1.6:           "<<Sta_OSevt_Eta[2]<<endl;

 double qmisid_eta_1,qmisid_eta_2,qmisid_eta_3,qmisid_err_1,qmisid_err_2,qmisid_err_3;
 qmisid_eta_1 = (double)Sta_SSevt_Eta[0]/(2*((double)Sta_SSevt_Eta[0] + (double)Sta_OSevt_Eta[0]));
 qmisid_err_1 = 0.5 * sqrt(((double)Sta_SSevt_Eta[0]+1)*((double)Sta_OSevt_Eta[0]+1)/(((double)Sta_SSevt_Eta[0] + (double)Sta_OSevt_Eta[0]+2)*
                     ((double)Sta_SSevt_Eta[0] + (double)Sta_OSevt_Eta[0]+2)*((double)Sta_SSevt_Eta[0] + (double)Sta_OSevt_Eta[0]+3)));
 qmisid_eta_2 = ((double)Sta_SSevt_Eta[2]/((double)Sta_SSevt_Eta[2]+(double)Sta_OSevt_Eta[2]) - qmisid_eta_1) / (1-2*qmisid_eta_1);

 qmisid_eta_3 = ((double)Sta_SSevt_Eta[2]/((double)Sta_SSevt_Eta[2]+(double)Sta_OSevt_Eta[2]) - qmisid_eta_1) / (1-2*qmisid_eta_1);

 double transfer_eta_2,transfer_eta_3;
 transfer_eta_2 = 0.5 * sqrt(((double)Sta_SSevt_Eta[1]+1)*((double)Sta_OSevt_Eta[1]+1)/(((double)Sta_SSevt_Eta[1] + (double)Sta_OSevt_Eta[1]+2)*
                     ((double)Sta_SSevt_Eta[1] + (double)Sta_OSevt_Eta[1]+2)*((double)Sta_SSevt_Eta[1] + (double)Sta_OSevt_Eta[1]+3)));
 transfer_eta_3 = 0.5 * sqrt(((double)Sta_SSevt_Eta[2]+1)*((double)Sta_OSevt_Eta[2]+1)/(((double)Sta_SSevt_Eta[2] + (double)Sta_OSevt_Eta[2]+2)*
                     ((double)Sta_SSevt_Eta[2] + (double)Sta_OSevt_Eta[2]+2)*((double)Sta_SSevt_Eta[2] + (double)Sta_OSevt_Eta[2]+3)));
 qmisid_err_2 = sqrt(transfer_eta_2*transfer_eta_2 * 1.0/((1-2*qmisid_eta_1)*(1-2*qmisid_eta_1)) +
                           qmisid_err_1*qmisid_err_1 * pow((2*(transfer_eta_2 - qmisid_eta_1)/((1-2*qmisid_eta_1)*(1-2*qmisid_eta_1)) - qmisid_eta_1/(1-2*qmisid_eta_1)), 2));

 qmisid_err_3 = sqrt(transfer_eta_3*transfer_eta_3 * 1.0/((1-2*qmisid_eta_1)*(1-2*qmisid_eta_1)) +
                           qmisid_err_1*qmisid_err_1 * pow((2*(transfer_eta_3 - qmisid_eta_1)/((1-2*qmisid_eta_1)*(1-2*qmisid_eta_1)) - qmisid_eta_1/(1-2*qmisid_eta_1)), 2));

 cout<<"q mis-id in |eta|<1.0:        "<<qmisid_eta_1<<" +/- "<<qmisid_err_1<<endl;
 cout<<"q mis-id in 1.0<|eta|<1.6:    "<<qmisid_eta_2<<" +/- "<<qmisid_err_2<<endl;
 cout<<"q mis-id in |eta|>1.6:        "<<qmisid_eta_3<<" +/- "<<qmisid_err_3<<endl;
*/
cout<<endl;
cout<<endl;
}

void Analysis::Initial(const char* RootName, int RootNumber)
{
//link to tree
  TTree *tree=NULL;
  TFile *file = (TFile *)gROOT->GetListOfFiles()->FindObject(RootName);
   if(!file)
    {file = new TFile(RootName);}
  tree = (TTree *)gDirectory->Get("t");
  if(tree == NULL){
    cout<<"No Such Tree!!!"<<endl;
  }else  cout<<"Tree succeeded taking out"<<endl;
  fChain = tree;
  //output: got a new rootfile
  cout<<"**Runing: Starting Rootfile "<<RootNumber<<endl;
  isRead=true;
  double filenumber;
  filenumber=fChain->GetEntries();
  if(filenumber==0)
    isRead=false;

//*****set address****//
// fChain->SetBranchAddress("zdata",&mc_evtwt);
    //mc global
    //
    fChain->SetBranchAddress("run",&run);
    fChain->SetBranchAddress("lumi",&lumi);
    fChain->SetBranchAddress("evt",&evt);
    fChain->SetBranchAddress("isData",&isData);
    fChain->SetBranchAddress("evt_scale1fb",&evt_scale1fb);
    fChain->SetBranchAddress("genps_weight",&genps_weight);
    fChain->SetBranchAddress("xsec_br",&xsec_br);
//    fChain->SetBranchAddress("CMS4path",&CMS4path);
    fChain->SetBranchAddress("CMS4index",&CMS4index);
    fChain->SetBranchAddress("weight_fr_r1_f1",&weight_fr_r1_f1);
    fChain->SetBranchAddress("weight_fr_r1_f2",&weight_fr_r1_f2);
    fChain->SetBranchAddress("weight_fr_r1_f0p5",&weight_fr_r1_f0p5);
    fChain->SetBranchAddress("weight_fr_r2_f1",&weight_fr_r2_f1);
    fChain->SetBranchAddress("weight_fr_r2_f2",&weight_fr_r2_f2);
    fChain->SetBranchAddress("weight_fr_r2_f0p5",&weight_fr_r2_f0p5);
    fChain->SetBranchAddress("weight_fr_r0p5_f1",&weight_fr_r0p5_f1);
    fChain->SetBranchAddress("weight_fr_r0p5_f2",&weight_fr_r0p5_f2);
    fChain->SetBranchAddress("weight_fr_r0p5_f0p5",&weight_fr_r0p5_f0p5);
    fChain->SetBranchAddress("weight_pdf_up",&weight_pdf_up);
    fChain->SetBranchAddress("weight_pdf_down",&weight_pdf_down);
    fChain->SetBranchAddress("weight_alphas_down",&weight_alphas_down);
    fChain->SetBranchAddress("weight_alphas_up",&weight_alphas_up);
    fChain->SetBranchAddress("HLT_DoubleMu",&HLT_DoubleMu);
    fChain->SetBranchAddress("HLT_DoubleEl",&HLT_DoubleEl);
    fChain->SetBranchAddress("HLT_MuEG",&HLT_MuEG);
    fChain->SetBranchAddress("pass_duplicate_ee_em_mm",&pass_duplicate_ee_em_mm);
    fChain->SetBranchAddress("pass_duplicate_mm_em_ee",&pass_duplicate_mm_em_ee);
    fChain->SetBranchAddress("gen_ht",&gen_ht);
    fChain->SetBranchAddress("firstgoodvertex",&firstgoodvertex);
    fChain->SetBranchAddress("nvtx",&nvtx);
    fChain->SetBranchAddress("nTrueInt",&nTrueInt);
//    fChain->SetBranchAddress("lep_p4_",&lep_p4_);
//    fChain->SetBranchAddress("lep_p4_fCoordinates_fX",&lep_p4_fCoordinates_fX);
//    fChain->SetBranchAddress("lep_p4_fCoordinates_fY",&lep_p4_fCoordinates_fY);
//    fChain->SetBranchAddress("lep_p4_fCoordinates_fZ",&lep_p4_fCoordinates_fZ);
//    fChain->SetBranchAddress("lep_p4_fCoordinates_fT",&lep_p4_fCoordinates_fT);
    fChain->SetBranchAddress("lep_pt",&lep_pt);
    fChain->SetBranchAddress("lep_eta",&lep_eta);
    fChain->SetBranchAddress("lep_phi",&lep_phi);
    fChain->SetBranchAddress("lep_energy",&lep_energy);
    fChain->SetBranchAddress("lep_mva",&lep_mva);
    fChain->SetBranchAddress("lep_relIso03EA",&lep_relIso03EA);
    fChain->SetBranchAddress("lep_relIso03EAwLep",&lep_relIso03EAwLep);
    fChain->SetBranchAddress("lep_ip3d",&lep_ip3d);
    fChain->SetBranchAddress("lep_sip3d",&lep_sip3d);
    fChain->SetBranchAddress("lep_dxy",&lep_dxy);
    fChain->SetBranchAddress("lep_dz",&lep_dz);
    fChain->SetBranchAddress("lep_mc_id",&lep_mc_id);
    fChain->SetBranchAddress("lep_motherIdv2",&lep_motherIdv2);
    fChain->SetBranchAddress("lep_idx",&lep_idx);
    fChain->SetBranchAddress("lep_id",&lep_id);
    fChain->SetBranchAddress("lep_isTightPOG",&lep_isTightPOG);
    fChain->SetBranchAddress("lep_isMediumPOG",&lep_isMediumPOG);

//    fChain->SetBranchAddress("fCoordinates_fX",&fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates_fY",&fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates_fZ",&fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates_fT",&fCoordinates_fT);
    fChain->SetBranchAddress("met_pt",&met_pt);
    fChain->SetBranchAddress("met_phi",&met_phi);
    fChain->SetBranchAddress("met_up_pt",&met_up_pt);
    fChain->SetBranchAddress("met_up_phi",&met_up_phi);
    fChain->SetBranchAddress("met_dn_pt",&met_dn_pt);
    fChain->SetBranchAddress("met_dn_phi",&met_dn_phi);
    fChain->SetBranchAddress("met_gen_pt",&met_gen_pt);
    fChain->SetBranchAddress("met_gen_phi",&met_gen_phi);
//    fChain->SetBranchAddress("jets_p4_fCoordinates_fX",&jets_p4_fCoordinates_fX);
//    fChain->SetBranchAddress("jets_p4_fCoordinates_fY",&jets_p4_fCoordinates_fY);
//    fChain->SetBranchAddress("jets_p4_fCoordinates_fZ",&jets_p4_fCoordinates_fZ);
//    fChain->SetBranchAddress("jets_p4_fCoordinates_fT",&jets_p4_fCoordinates_fT);
    fChain->SetBranchAddress("jets_pt",&jets_pt);
    fChain->SetBranchAddress("jets_eta",&jets_eta);
    fChain->SetBranchAddress("jets_phi",&jets_phi);
    fChain->SetBranchAddress("jets_mass",&jets_mass);
    fChain->SetBranchAddress("nj",&nj);
    fChain->SetBranchAddress("nb",&nb);
    fChain->SetBranchAddress("nbmed",&nbmed);
    fChain->SetBranchAddress("ht",&ht);
    fChain->SetBranchAddress("weight_btagsf",&weight_btagsf);
    fChain->SetBranchAddress("weight_btagsf_heavy_DN",&weight_btagsf_heavy_DN);
    fChain->SetBranchAddress("weight_btagsf_heavy_UP",&weight_btagsf_heavy_UP);
    fChain->SetBranchAddress("weight_btagsf_light_DN",&weight_btagsf_light_DN);
    fChain->SetBranchAddress("weight_btagsf_light_UP",&weight_btagsf_light_UP);

/*    fChain->SetBranchAddress("",&);
    fChain->SetBranchAddress("",&);
*/



  
}

#endif //#ifdef zeeAna_cxx

