//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 22 18:20:33 2019 by ROOT version 6.08/07
// from TTree t_www/Signal Events
// found on file: results.root
//////////////////////////////////////////////////////////

#ifndef t_www_h
#define t_www_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TString.h"
#include "vector"
#include "Math/GenVector/PxPyPzE4D.h"
#include "vector"
#include "vector"
#include "Math/GenVector/LorentzVector.h"

class t_www {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxlep_p4 = 5;
   const Int_t kMaxjets_p4 = 11;
   const Int_t kMaxjets_up_p4 = 12;
   const Int_t kMaxjets_dn_p4 = 11;
   const Int_t kMaxjets_jer_p4 = 22;
   const Int_t kMaxjets_jerup_p4 = 23;
   const Int_t kMaxjets_jerdn_p4 = 13;
   const Int_t kMaxjets30_p4 = 6;
   const Int_t kMaxjets30_up_p4 = 6;
   const Int_t kMaxjets30_dn_p4 = 6;
   const Int_t kMaxjets30_jer_p4 = 6;
   const Int_t kMaxjets30_jerup_p4 = 6;
   const Int_t kMaxjets30_jerdn_p4 = 6;
   const Int_t kMaxak8jets_p4 = 3;
   const Int_t kMaxgenPart_p4 = 220;
   const Int_t kMaxdecay_p4 = 2;
   const Int_t kMaxlepton_p4 = 2;
   const Int_t kMaxquark_p4 = 4;
   const Int_t kMaxboosted0_decay_p4 = 2;
   const Int_t kMaxboosted0_lepton_p4 = 2;
   const Int_t kMaxboosted0_quark_p4 = 4;
   const Int_t kMaxboosted250_decay_p4 = 2;
   const Int_t kMaxboosted250_lepton_p4 = 2;
   const Int_t kMaxboosted250_quark_p4 = 4;
   const Int_t kMaxboosted500_decay_p4 = 2;
   const Int_t kMaxboosted500_lepton_p4 = 2;
   const Int_t kMaxboosted500_quark_p4 = 4;
   const Int_t kMaxboosted1000_decay_p4 = 2;
   const Int_t kMaxboosted1000_lepton_p4 = 2;
   const Int_t kMaxboosted1000_quark_p4 = 4;
   const Int_t kMaxboosted1500_decay_p4 = 2;
   const Int_t kMaxboosted1500_lepton_p4 = 2;
   const Int_t kMaxboosted1500_quark_p4 = 4;
   const Int_t kMaxw_p4 = 3;
   const Int_t kMaxl_p4 = 5;
   const Int_t kMaxq_p4 = 4;

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   ULong64_t       evt;
   Int_t           isData;
   Float_t         evt_scale1fb;
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
   Float_t         weight_isr;
   Float_t         weight_isr_up;
   Float_t         weight_isr_down;
   Int_t           HLT_DoubleMu;
   Int_t           HLT_DoubleEl;
   Int_t           HLT_DoubleEl_DZ;
   Int_t           HLT_DoubleEl_DZ_2;
   Int_t           HLT_MuEG;
   Int_t           HLT_SingleEl8;
   Int_t           HLT_SingleEl17;
   Int_t           HLT_SingleIsoEl8;
   Int_t           HLT_SingleIsoEl17;
   Int_t           HLT_SingleIsoEl23;
   Int_t           HLT_SingleIsoMu8;
   Int_t           HLT_SingleIsoMu17;
   Int_t           HLT_PFMET140_PFMHT140_IDTight;
   Int_t           mc_HLT_DoubleMu;
   Int_t           mc_HLT_DoubleEl;
   Int_t           mc_HLT_DoubleEl_DZ;
   Int_t           mc_HLT_DoubleEl_DZ_2;
   Int_t           mc_HLT_MuEG;
   Int_t           mc_HLT_SingleEl8;
   Int_t           mc_HLT_SingleEl17;
   Int_t           mc_HLT_SingleIsoEl8;
   Int_t           mc_HLT_SingleIsoEl17;
   Int_t           mc_HLT_SingleIsoEl23;
   Int_t           mc_HLT_SingleIsoMu8;
   Int_t           mc_HLT_SingleIsoMu17;
   Int_t           mc_HLT_PFMET140_PFMHT140_IDTight;
   Int_t           pass_duplicate_ee_em_mm;
   Int_t           pass_duplicate_mm_em_ee;
   Int_t           is2016;
   Int_t           is2017;
   Int_t           HLT_MuEG_2016;
   Int_t           mc_HLT_MuEG_2016;
   Int_t           pass_duplicate_ee_em2016_mm;
   Int_t           pass_duplicate_mm_em2016_ee;
   Int_t           passTrigger;
   Int_t           lep_p4_;
   Float_t         lep_p4_fCoordinates_fX[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fY[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fZ[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fT[kMaxlep_p4];   //[lep_p4_]
   vector<float>   *lep_pt;
   vector<float>   *lep_eta;
   vector<float>   *lep_phi;
   vector<float>   *lep_coneCorrPt;
   vector<float>   *lep_ip3d;
   vector<float>   *lep_ip3derr;
   vector<int>     *lep_isTriggerSafe_v1;
   vector<int>     *lep_lostHits;
   vector<int>     *lep_convVeto;
   vector<int>     *lep_motherIdSS;
   vector<int>     *lep_pass_VVV_cutbased_3l_fo;
   vector<int>     *lep_pass_VVV_cutbased_3l_tight;
   vector<int>     *lep_pass_VVV_cutbased_fo;
   vector<int>     *lep_pass_VVV_cutbased_tight;
   vector<int>     *lep_pass_VVV_cutbased_veto;
   vector<int>     *lep_pass_VVV_cutbased_fo_noiso;
   vector<int>     *lep_pass_VVV_cutbased_tight_noiso;
   vector<int>     *lep_pass_VVV_cutbased_veto_noiso;
   vector<int>     *lep_pass_POG_veto;
   vector<int>     *lep_pass_POG_loose;
   vector<int>     *lep_pass_POG_medium;
   vector<int>     *lep_pass_POG_tight;
   vector<int>     *lep_pdgId;
   vector<float>   *lep_dxy;
   vector<float>   *lep_dz;
   vector<float>   *lep_pterr;
   vector<float>   *lep_relIso04DB;
   vector<float>   *lep_relIso03EA;
   vector<float>   *lep_relIso03EALep;
   vector<float>   *lep_relIso03EAv2;
   vector<float>   *lep_relIso04EAv2;
   vector<float>   *lep_relIso03EAv2Lep;
   vector<int>     *lep_tightCharge;
   vector<float>   *lep_trk_pt;
   vector<int>     *lep_charge;
   vector<float>   *lep_etaSC;
   vector<float>   *lep_MVA;
   vector<int>     *lep_isMediumPOG;
   vector<int>     *lep_isTightPOG;
   vector<int>     *lep_isFromW;
   vector<int>     *lep_isFromZ;
   vector<int>     *lep_isFromB;
   vector<int>     *lep_isFromC;
   vector<int>     *lep_isFromL;
   vector<int>     *lep_isFromLF;
   vector<int>     *lep_genPart_index;
   vector<float>   *lep_r9;
   vector<int>     *lep_nlayers;
   Float_t         el_pt;
   Float_t         el_eta;
   Float_t         el_phi;
   Float_t         el_relIso03EA;
   Float_t         el_relIso03EALep;
   Float_t         el_ip3d;
   Float_t         mu_pt;
   Float_t         mu_eta;
   Float_t         mu_phi;
   Float_t         mu_relIso04DB;
   Float_t         mu_relIso03EA;
   Float_t         mu_relIso03EALep;
   Float_t         mu_ip3d;
   Float_t         lbnt_pt;
   Float_t         lbnt_coneCorrPt;
   Float_t         lbnt_abseta;
   Float_t         lbnt_pdgId;
   Float_t         lbnt_el_pt;
   Float_t         lbnt_el_coneCorrPt;
   Float_t         lbnt_el_abseta;
   Float_t         lbnt_mu_pt;
   Float_t         lbnt_mu_coneCorrPt;
   Float_t         lbnt_mu_abseta;
   Int_t           jets_p4_;
   Float_t         jets_p4_fCoordinates_fX[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fY[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fZ[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fT[kMaxjets_p4];   //[jets_p4_]
   Int_t           jets_up_p4_;
   Float_t         jets_up_p4_fCoordinates_fX[kMaxjets_up_p4];   //[jets_up_p4_]
   Float_t         jets_up_p4_fCoordinates_fY[kMaxjets_up_p4];   //[jets_up_p4_]
   Float_t         jets_up_p4_fCoordinates_fZ[kMaxjets_up_p4];   //[jets_up_p4_]
   Float_t         jets_up_p4_fCoordinates_fT[kMaxjets_up_p4];   //[jets_up_p4_]
   Int_t           jets_dn_p4_;
   Float_t         jets_dn_p4_fCoordinates_fX[kMaxjets_dn_p4];   //[jets_dn_p4_]
   Float_t         jets_dn_p4_fCoordinates_fY[kMaxjets_dn_p4];   //[jets_dn_p4_]
   Float_t         jets_dn_p4_fCoordinates_fZ[kMaxjets_dn_p4];   //[jets_dn_p4_]
   Float_t         jets_dn_p4_fCoordinates_fT[kMaxjets_dn_p4];   //[jets_dn_p4_]
   vector<float>   *jets_csv;
   vector<float>   *jets_up_csv;
   vector<float>   *jets_dn_csv;
   vector<float>   *jets_jer_csv;
   vector<float>   *jets_jerup_csv;
   vector<float>   *jets_jerdn_csv;
   Int_t           jets_jer_p4_;
   Float_t         jets_jer_p4_fCoordinates_fX[kMaxjets_jer_p4];   //[jets_jer_p4_]
   Float_t         jets_jer_p4_fCoordinates_fY[kMaxjets_jer_p4];   //[jets_jer_p4_]
   Float_t         jets_jer_p4_fCoordinates_fZ[kMaxjets_jer_p4];   //[jets_jer_p4_]
   Float_t         jets_jer_p4_fCoordinates_fT[kMaxjets_jer_p4];   //[jets_jer_p4_]
   Int_t           jets_jerup_p4_;
   Float_t         jets_jerup_p4_fCoordinates_fX[kMaxjets_jerup_p4];   //[jets_jerup_p4_]
   Float_t         jets_jerup_p4_fCoordinates_fY[kMaxjets_jerup_p4];   //[jets_jerup_p4_]
   Float_t         jets_jerup_p4_fCoordinates_fZ[kMaxjets_jerup_p4];   //[jets_jerup_p4_]
   Float_t         jets_jerup_p4_fCoordinates_fT[kMaxjets_jerup_p4];   //[jets_jerup_p4_]
   Int_t           jets_jerdn_p4_;
   Float_t         jets_jerdn_p4_fCoordinates_fX[kMaxjets_jerdn_p4];   //[jets_jerdn_p4_]
   Float_t         jets_jerdn_p4_fCoordinates_fY[kMaxjets_jerdn_p4];   //[jets_jerdn_p4_]
   Float_t         jets_jerdn_p4_fCoordinates_fZ[kMaxjets_jerdn_p4];   //[jets_jerdn_p4_]
   Float_t         jets_jerdn_p4_fCoordinates_fT[kMaxjets_jerdn_p4];   //[jets_jerdn_p4_]
   Int_t           jets30_p4_;
   Float_t         jets30_p4_fCoordinates_fX[kMaxjets30_p4];   //[jets30_p4_]
   Float_t         jets30_p4_fCoordinates_fY[kMaxjets30_p4];   //[jets30_p4_]
   Float_t         jets30_p4_fCoordinates_fZ[kMaxjets30_p4];   //[jets30_p4_]
   Float_t         jets30_p4_fCoordinates_fT[kMaxjets30_p4];   //[jets30_p4_]
   Int_t           jets30_up_p4_;
   Float_t         jets30_up_p4_fCoordinates_fX[kMaxjets30_up_p4];   //[jets30_up_p4_]
   Float_t         jets30_up_p4_fCoordinates_fY[kMaxjets30_up_p4];   //[jets30_up_p4_]
   Float_t         jets30_up_p4_fCoordinates_fZ[kMaxjets30_up_p4];   //[jets30_up_p4_]
   Float_t         jets30_up_p4_fCoordinates_fT[kMaxjets30_up_p4];   //[jets30_up_p4_]
   Int_t           jets30_dn_p4_;
   Float_t         jets30_dn_p4_fCoordinates_fX[kMaxjets30_dn_p4];   //[jets30_dn_p4_]
   Float_t         jets30_dn_p4_fCoordinates_fY[kMaxjets30_dn_p4];   //[jets30_dn_p4_]
   Float_t         jets30_dn_p4_fCoordinates_fZ[kMaxjets30_dn_p4];   //[jets30_dn_p4_]
   Float_t         jets30_dn_p4_fCoordinates_fT[kMaxjets30_dn_p4];   //[jets30_dn_p4_]
   Int_t           jets30_jer_p4_;
   Float_t         jets30_jer_p4_fCoordinates_fX[kMaxjets30_jer_p4];   //[jets30_jer_p4_]
   Float_t         jets30_jer_p4_fCoordinates_fY[kMaxjets30_jer_p4];   //[jets30_jer_p4_]
   Float_t         jets30_jer_p4_fCoordinates_fZ[kMaxjets30_jer_p4];   //[jets30_jer_p4_]
   Float_t         jets30_jer_p4_fCoordinates_fT[kMaxjets30_jer_p4];   //[jets30_jer_p4_]
   Int_t           jets30_jerup_p4_;
   Float_t         jets30_jerup_p4_fCoordinates_fX[kMaxjets30_jerup_p4];   //[jets30_jerup_p4_]
   Float_t         jets30_jerup_p4_fCoordinates_fY[kMaxjets30_jerup_p4];   //[jets30_jerup_p4_]
   Float_t         jets30_jerup_p4_fCoordinates_fZ[kMaxjets30_jerup_p4];   //[jets30_jerup_p4_]
   Float_t         jets30_jerup_p4_fCoordinates_fT[kMaxjets30_jerup_p4];   //[jets30_jerup_p4_]
   Int_t           jets30_jerdn_p4_;
   Float_t         jets30_jerdn_p4_fCoordinates_fX[kMaxjets30_jerdn_p4];   //[jets30_jerdn_p4_]
   Float_t         jets30_jerdn_p4_fCoordinates_fY[kMaxjets30_jerdn_p4];   //[jets30_jerdn_p4_]
   Float_t         jets30_jerdn_p4_fCoordinates_fZ[kMaxjets30_jerdn_p4];   //[jets30_jerdn_p4_]
   Float_t         jets30_jerdn_p4_fCoordinates_fT[kMaxjets30_jerdn_p4];   //[jets30_jerdn_p4_]
   Int_t           ak8jets_p4_;
   Float_t         ak8jets_p4_fCoordinates_fX[kMaxak8jets_p4];   //[ak8jets_p4_]
   Float_t         ak8jets_p4_fCoordinates_fY[kMaxak8jets_p4];   //[ak8jets_p4_]
   Float_t         ak8jets_p4_fCoordinates_fZ[kMaxak8jets_p4];   //[ak8jets_p4_]
   Float_t         ak8jets_p4_fCoordinates_fT[kMaxak8jets_p4];   //[ak8jets_p4_]
   vector<float>   *ak8jets_softdropMass;
   vector<float>   *ak8jets_prunedMass;
   vector<float>   *ak8jets_trimmedMass;
   vector<float>   *ak8jets_mass;
   vector<float>   *ak8jets_nJettinessTau1;
   vector<float>   *ak8jets_nJettinessTau2;
   vector<float>   *ak8jets_softdropPuppiSubjet1;
   vector<float>   *ak8jets_softdropPuppiSubjet2;
   vector<float>   *ak8jets_puppi_softdropMass;
   vector<float>   *ak8jets_puppi_nJettinessTau1;
   vector<float>   *ak8jets_puppi_nJettinessTau2;
   vector<float>   *ak8jets_puppi_eta;
   vector<float>   *ak8jets_puppi_phi;
   vector<float>   *ak8jets_puppi_pt;
   vector<float>   *ak8jets_puppi_mass;
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         met_up_pt;
   Float_t         met_up_phi;
   Float_t         met_dn_pt;
   Float_t         met_dn_phi;
   Float_t         met_gen_pt;
   Float_t         met_gen_phi;
   Float_t         met_jer_pt;
   Float_t         met_jerup_pt;
   Float_t         met_jerdn_pt;
   Float_t         met_jer_phi;
   Float_t         met_jerup_phi;
   Float_t         met_jerdn_phi;
   Int_t           firstgoodvertex;
   Int_t           nTrueInt;
   Int_t           nVert;
   Int_t           nisoTrack_mt2_cleaned_VVV_cutbased_veto;
   Float_t         weight_btagsf;
   Float_t         weight_btagsf_heavy_DN;
   Float_t         weight_btagsf_heavy_UP;
   Float_t         weight_btagsf_light_DN;
   Float_t         weight_btagsf_light_UP;
   Float_t         gen_ht;
   Int_t           genPart_p4_;
   Float_t         genPart_p4_fCoordinates_fX[kMaxgenPart_p4];   //[genPart_p4_]
   Float_t         genPart_p4_fCoordinates_fY[kMaxgenPart_p4];   //[genPart_p4_]
   Float_t         genPart_p4_fCoordinates_fZ[kMaxgenPart_p4];   //[genPart_p4_]
   Float_t         genPart_p4_fCoordinates_fT[kMaxgenPart_p4];   //[genPart_p4_]
   vector<int>     *genPart_motherId;
   vector<int>     *genPart_pdgId;
   vector<int>     *genPart_charge;
   vector<int>     *genPart_status;
   Int_t           ngenLep;
   Int_t           ngenLepFromTau;
   Int_t           ngenLepFromBoson;
   Int_t           Flag_AllEventFilters;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_HBHEIsoNoiseFilter;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_badChargedCandidateFilter;
   Int_t           Flag_badMuonFilter;
   Int_t           Flag_badMuonFilterv2;
   Int_t           Flag_badChargedCandidateFilterv2;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_ecalBadCalibFilter;
   Int_t           Flag_globalTightHalo2016;
   Int_t           Flag_goodVertices;
   Int_t           Flag_ecalLaserCorrFilter;
   Int_t           Flag_hcalLaserEventFilter;
   Int_t           Flag_trackingFailureFilter;
   Int_t           Flag_CSCTightHaloFilter;
   Int_t           Flag_CSCTightHalo2015Filter;
   Int_t           Flag_badMuons;
   Int_t           Flag_duplicateMuons;
   Int_t           Flag_noBadMuons;
   Int_t           fastsimfilt;
   Int_t           nVlep;
   Int_t           nTlep;
   Int_t           nTlepSS;
   Int_t           nLlep;
   Int_t           nLlep3L;
   Int_t           nTlep3L;
   Int_t           nSFOS;
   Int_t           nSFOSinZ;
   Int_t           nj;
   Int_t           nj_up;
   Int_t           nj_dn;
   Int_t           nj_jer;
   Int_t           nj_jerup;
   Int_t           nj_jerdn;
   Int_t           nj30;
   Int_t           nj30_up;
   Int_t           nj30_dn;
   Int_t           nj30_jer;
   Int_t           nj30_jerup;
   Int_t           nj30_jerdn;
   Int_t           nb;
   Int_t           nb_up;
   Int_t           nb_dn;
   Int_t           nb_jer;
   Int_t           nb_jerup;
   Int_t           nb_jerdn;
   Float_t         Ml0j0;
   Float_t         Ml0j0_up;
   Float_t         Ml0j0_dn;
   Float_t         Ml0j0_jer;
   Float_t         Ml0j0_jerup;
   Float_t         Ml0j0_jerdn;
   Float_t         Ml0j1;
   Float_t         Ml0j1_up;
   Float_t         Ml0j1_dn;
   Float_t         Ml0j1_jer;
   Float_t         Ml0j1_jerup;
   Float_t         Ml0j1_jerdn;
   Float_t         Ml1j0;
   Float_t         Ml1j0_up;
   Float_t         Ml1j0_dn;
   Float_t         Ml1j0_jer;
   Float_t         Ml1j0_jerup;
   Float_t         Ml1j0_jerdn;
   Float_t         Ml1j1;
   Float_t         Ml1j1_up;
   Float_t         Ml1j1_dn;
   Float_t         Ml1j1_jer;
   Float_t         Ml1j1_jerup;
   Float_t         Ml1j1_jerdn;
   Float_t         MinMlj;
   Float_t         MinMlj_up;
   Float_t         MinMlj_dn;
   Float_t         MinMlj_jer;
   Float_t         MinMlj_jerup;
   Float_t         MinMlj_jerdn;
   Float_t         SumMinMlj01;
   Float_t         SumMinMlj01_up;
   Float_t         SumMinMlj01_dn;
   Float_t         SumMinMlj01_jer;
   Float_t         SumMinMlj01_jerup;
   Float_t         SumMinMlj01_jerdn;
   Float_t         MaxMlj;
   Float_t         MaxMlj_up;
   Float_t         MaxMlj_dn;
   Float_t         MaxMlj_jer;
   Float_t         MaxMlj_jerup;
   Float_t         MaxMlj_jerdn;
   Float_t         SumMlj;
   Float_t         SumMlj_up;
   Float_t         SumMlj_dn;
   Float_t         SumMlj_jer;
   Float_t         SumMlj_jerup;
   Float_t         SumMlj_jerdn;
   Float_t         Ml0jj;
   Float_t         Ml0jj_up;
   Float_t         Ml0jj_dn;
   Float_t         Ml0jj_jer;
   Float_t         Ml0jj_jerup;
   Float_t         Ml0jj_jerdn;
   Float_t         Ml1jj;
   Float_t         Ml1jj_up;
   Float_t         Ml1jj_dn;
   Float_t         Ml1jj_jer;
   Float_t         Ml1jj_jerup;
   Float_t         Ml1jj_jerdn;
   Float_t         MinMljj;
   Float_t         MinMljj_up;
   Float_t         MinMljj_dn;
   Float_t         MinMljj_jer;
   Float_t         MinMljj_jerup;
   Float_t         MinMljj_jerdn;
   Float_t         MaxMljj;
   Float_t         MaxMljj_up;
   Float_t         MaxMljj_dn;
   Float_t         MaxMljj_jer;
   Float_t         MaxMljj_jerup;
   Float_t         MaxMljj_jerdn;
   Float_t         SumMljj;
   Float_t         SumMljj_up;
   Float_t         SumMljj_dn;
   Float_t         SumMljj_jer;
   Float_t         SumMljj_jerup;
   Float_t         SumMljj_jerdn;
   Float_t         Mjj;
   Float_t         Mjj_up;
   Float_t         Mjj_dn;
   Float_t         Mjj_jer;
   Float_t         Mjj_jerup;
   Float_t         Mjj_jerdn;
   Float_t         DRjj;
   Float_t         DRjj_up;
   Float_t         DRjj_dn;
   Float_t         DRjj_jer;
   Float_t         DRjj_jerup;
   Float_t         DRjj_jerdn;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_up;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_dn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_jer;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_jerup;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_jerdn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_up;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_dn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_jer;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_jerup;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_jerdn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
   Float_t         MjjDR1;
   Float_t         MjjDR1_up;
   Float_t         MjjDR1_dn;
   Float_t         MjjDR1_jer;
   Float_t         MjjDR1_jerup;
   Float_t         MjjDR1_jerdn;
   Float_t         DRjjDR1;
   Float_t         DRjjDR1_up;
   Float_t         DRjjDR1_dn;
   Float_t         DRjjDR1_jer;
   Float_t         DRjjDR1_jerup;
   Float_t         DRjjDR1_jerdn;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1_up;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1_dn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1_jer;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1_jerup;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet0_wtag_p4_DR1_jerdn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1_up;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1_dn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1_jer;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1_jerup;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *jet1_wtag_p4_DR1_jerdn;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
   Float_t         MjjVBF;
   Float_t         MjjVBF_up;
   Float_t         MjjVBF_dn;
   Float_t         MjjVBF_jer;
   Float_t         MjjVBF_jerup;
   Float_t         MjjVBF_jerdn;
   Float_t         DetajjVBF;
   Float_t         DetajjVBF_up;
   Float_t         DetajjVBF_dn;
   Float_t         DetajjVBF_jer;
   Float_t         DetajjVBF_jerup;
   Float_t         DetajjVBF_jerdn;
   Float_t         MjjL;
   Float_t         MjjL_up;
   Float_t         MjjL_dn;
   Float_t         MjjL_jer;
   Float_t         MjjL_jerup;
   Float_t         MjjL_jerdn;
   Float_t         DetajjL;
   Float_t         DetajjL_up;
   Float_t         DetajjL_dn;
   Float_t         DetajjL_jer;
   Float_t         DetajjL_jerup;
   Float_t         DetajjL_jerdn;
   Float_t         MllSS;
   Float_t         MeeSS;
   Float_t         Mll3L;
   Float_t         Mee3L;
   Float_t         Mll3L1;
   Float_t         M3l;
   Float_t         Pt3l;
   Float_t         M01;
   Float_t         M02;
   Float_t         M12;
   Int_t           isSFOS01;
   Int_t           isSFOS02;
   Int_t           isSFOS12;
   Float_t         DPhi3lMET;
   Float_t         DPhi3lMET_up;
   Float_t         DPhi3lMET_dn;
   Float_t         DPhi3lMET_jer;
   Float_t         DPhi3lMET_jerup;
   Float_t         DPhi3lMET_jerdn;
   Float_t         DPhi3lMET_gen;
   Float_t         MTmax;
   Float_t         MTmax_up;
   Float_t         MTmax_dn;
   Float_t         MTmax_jer;
   Float_t         MTmax_jerup;
   Float_t         MTmax_jerdn;
   Float_t         MTmax_gen;
   Float_t         MTmin;
   Float_t         MTmin_up;
   Float_t         MTmin_dn;
   Float_t         MTmin_jer;
   Float_t         MTmin_jerup;
   Float_t         MTmin_jerdn;
   Float_t         MTmin_gen;
   Float_t         MT3rd;
   Float_t         MT3rd_up;
   Float_t         MT3rd_dn;
   Float_t         MT3rd_jer;
   Float_t         MT3rd_jerup;
   Float_t         MT3rd_jerdn;
   Float_t         MT3rd_gen;
   Float_t         MTmax3L;
   Float_t         MTmax3L_up;
   Float_t         MTmax3L_dn;
   Float_t         MTmax3L_jer;
   Float_t         MTmax3L_jerup;
   Float_t         MTmax3L_jerdn;
   Float_t         MTmax3L_gen;
   Int_t           passSSee;
   Int_t           passSSem;
   Int_t           passSSmm;
   Int_t           lep_idx0_SS;
   Int_t           lep_idx1_SS;
   TString         *bkgtype;
   Int_t           vetophoton;
   Float_t         purewgt;
   Float_t         purewgt_up;
   Float_t         purewgt_dn;
   Float_t         ffwgt;
   Float_t         ffwgt_up;
   Float_t         ffwgt_dn;
   Float_t         ffwgt_el_up;
   Float_t         ffwgt_el_dn;
   Float_t         ffwgt_mu_up;
   Float_t         ffwgt_mu_dn;
   Float_t         ffwgt_closure_up;
   Float_t         ffwgt_closure_dn;
   Float_t         ffwgt_closure_el_up;
   Float_t         ffwgt_closure_el_dn;
   Float_t         ffwgt_closure_mu_up;
   Float_t         ffwgt_closure_mu_dn;
   Float_t         ffwgt_full_up;
   Float_t         ffwgt_full_dn;
   Float_t         ffwgtqcd;
   Float_t         ffwgtqcd_up;
   Float_t         ffwgtqcd_dn;
   Float_t         lepsf;
   Float_t         lepsf_up;
   Float_t         lepsf_dn;
   Float_t         trigeff;
   Float_t         trigeff_up;
   Float_t         trigeff_dn;
   Float_t         trigsf;
   Float_t         trigsf_up;
   Float_t         trigsf_dn;
   Float_t         musmear_sf;
   Int_t           higgsdecay;
   Int_t           nlep;
   Int_t           nquark;
   Int_t           wa_id;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *higgs_p4;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
   Int_t           decay_p4_;
   Float_t         decay_p4_fCoordinates_fX[kMaxdecay_p4];   //[decay_p4_]
   Float_t         decay_p4_fCoordinates_fY[kMaxdecay_p4];   //[decay_p4_]
   Float_t         decay_p4_fCoordinates_fZ[kMaxdecay_p4];   //[decay_p4_]
   Float_t         decay_p4_fCoordinates_fT[kMaxdecay_p4];   //[decay_p4_]
   vector<int>     *decay_id;
   vector<int>     *decay_isstar;
   Int_t           lepton_p4_;
   Float_t         lepton_p4_fCoordinates_fX[kMaxlepton_p4];   //[lepton_p4_]
   Float_t         lepton_p4_fCoordinates_fY[kMaxlepton_p4];   //[lepton_p4_]
   Float_t         lepton_p4_fCoordinates_fZ[kMaxlepton_p4];   //[lepton_p4_]
   Float_t         lepton_p4_fCoordinates_fT[kMaxlepton_p4];   //[lepton_p4_]
   vector<int>     *lepton_id;
   vector<int>     *lepton_isstar;
   Int_t           quark_p4_;
   Float_t         quark_p4_fCoordinates_fX[kMaxquark_p4];   //[quark_p4_]
   Float_t         quark_p4_fCoordinates_fY[kMaxquark_p4];   //[quark_p4_]
   Float_t         quark_p4_fCoordinates_fZ[kMaxquark_p4];   //[quark_p4_]
   Float_t         quark_p4_fCoordinates_fT[kMaxquark_p4];   //[quark_p4_]
   vector<int>     *quark_id;
   vector<int>     *quark_isstar;
   Float_t         lqq_max_dr;
   Float_t         lq0_dr;
   Float_t         lq1_dr;
   Float_t         qq_dr;
   Int_t           boosted0_decay_p4_;
   Float_t         boosted0_decay_p4_fCoordinates_fX[kMaxboosted0_decay_p4];   //[boosted0_decay_p4_]
   Float_t         boosted0_decay_p4_fCoordinates_fY[kMaxboosted0_decay_p4];   //[boosted0_decay_p4_]
   Float_t         boosted0_decay_p4_fCoordinates_fZ[kMaxboosted0_decay_p4];   //[boosted0_decay_p4_]
   Float_t         boosted0_decay_p4_fCoordinates_fT[kMaxboosted0_decay_p4];   //[boosted0_decay_p4_]
   vector<int>     *boosted0_decay_id;
   vector<int>     *boosted0_decay_isstar;
   vector<float>   *boosted0_decay_h_dr;
   vector<float>   *boosted0_decay_h_deta;
   vector<float>   *boosted0_decay_h_dphi;
   vector<float>   *boosted0_decay_h_deta_rotated;
   vector<float>   *boosted0_decay_h_dphi_rotated;
   Int_t           boosted0_lepton_p4_;
   Float_t         boosted0_lepton_p4_fCoordinates_fX[kMaxboosted0_lepton_p4];   //[boosted0_lepton_p4_]
   Float_t         boosted0_lepton_p4_fCoordinates_fY[kMaxboosted0_lepton_p4];   //[boosted0_lepton_p4_]
   Float_t         boosted0_lepton_p4_fCoordinates_fZ[kMaxboosted0_lepton_p4];   //[boosted0_lepton_p4_]
   Float_t         boosted0_lepton_p4_fCoordinates_fT[kMaxboosted0_lepton_p4];   //[boosted0_lepton_p4_]
   vector<int>     *boosted0_lepton_id;
   vector<int>     *boosted0_lepton_isstar;
   vector<float>   *boosted0_lepton_h_dr;
   vector<float>   *boosted0_lepton_h_deta;
   vector<float>   *boosted0_lepton_h_dphi;
   vector<float>   *boosted0_lepton_h_deta_rotated;
   vector<float>   *boosted0_lepton_h_dphi_rotated;
   Int_t           boosted0_quark_p4_;
   Float_t         boosted0_quark_p4_fCoordinates_fX[kMaxboosted0_quark_p4];   //[boosted0_quark_p4_]
   Float_t         boosted0_quark_p4_fCoordinates_fY[kMaxboosted0_quark_p4];   //[boosted0_quark_p4_]
   Float_t         boosted0_quark_p4_fCoordinates_fZ[kMaxboosted0_quark_p4];   //[boosted0_quark_p4_]
   Float_t         boosted0_quark_p4_fCoordinates_fT[kMaxboosted0_quark_p4];   //[boosted0_quark_p4_]
   vector<int>     *boosted0_quark_id;
   vector<int>     *boosted0_quark_isstar;
   vector<float>   *boosted0_quark_h_dr;
   vector<float>   *boosted0_quark_h_deta;
   vector<float>   *boosted0_quark_h_dphi;
   vector<float>   *boosted0_quark_h_deta_rotated;
   vector<float>   *boosted0_quark_h_dphi_rotated;
   Int_t           boosted250_decay_p4_;
   Float_t         boosted250_decay_p4_fCoordinates_fX[kMaxboosted250_decay_p4];   //[boosted250_decay_p4_]
   Float_t         boosted250_decay_p4_fCoordinates_fY[kMaxboosted250_decay_p4];   //[boosted250_decay_p4_]
   Float_t         boosted250_decay_p4_fCoordinates_fZ[kMaxboosted250_decay_p4];   //[boosted250_decay_p4_]
   Float_t         boosted250_decay_p4_fCoordinates_fT[kMaxboosted250_decay_p4];   //[boosted250_decay_p4_]
   vector<int>     *boosted250_decay_id;
   vector<int>     *boosted250_decay_isstar;
   vector<float>   *boosted250_decay_h_dr;
   vector<float>   *boosted250_decay_h_deta;
   vector<float>   *boosted250_decay_h_dphi;
   vector<float>   *boosted250_decay_h_deta_rotated;
   vector<float>   *boosted250_decay_h_dphi_rotated;
   Int_t           boosted250_lepton_p4_;
   Float_t         boosted250_lepton_p4_fCoordinates_fX[kMaxboosted250_lepton_p4];   //[boosted250_lepton_p4_]
   Float_t         boosted250_lepton_p4_fCoordinates_fY[kMaxboosted250_lepton_p4];   //[boosted250_lepton_p4_]
   Float_t         boosted250_lepton_p4_fCoordinates_fZ[kMaxboosted250_lepton_p4];   //[boosted250_lepton_p4_]
   Float_t         boosted250_lepton_p4_fCoordinates_fT[kMaxboosted250_lepton_p4];   //[boosted250_lepton_p4_]
   vector<int>     *boosted250_lepton_id;
   vector<int>     *boosted250_lepton_isstar;
   vector<float>   *boosted250_lepton_h_dr;
   vector<float>   *boosted250_lepton_h_deta;
   vector<float>   *boosted250_lepton_h_dphi;
   vector<float>   *boosted250_lepton_h_deta_rotated;
   vector<float>   *boosted250_lepton_h_dphi_rotated;
   Int_t           boosted250_quark_p4_;
   Float_t         boosted250_quark_p4_fCoordinates_fX[kMaxboosted250_quark_p4];   //[boosted250_quark_p4_]
   Float_t         boosted250_quark_p4_fCoordinates_fY[kMaxboosted250_quark_p4];   //[boosted250_quark_p4_]
   Float_t         boosted250_quark_p4_fCoordinates_fZ[kMaxboosted250_quark_p4];   //[boosted250_quark_p4_]
   Float_t         boosted250_quark_p4_fCoordinates_fT[kMaxboosted250_quark_p4];   //[boosted250_quark_p4_]
   vector<int>     *boosted250_quark_id;
   vector<int>     *boosted250_quark_isstar;
   vector<float>   *boosted250_quark_h_dr;
   vector<float>   *boosted250_quark_h_deta;
   vector<float>   *boosted250_quark_h_dphi;
   vector<float>   *boosted250_quark_h_deta_rotated;
   vector<float>   *boosted250_quark_h_dphi_rotated;
   Float_t         boosted250_lqq_max_dr;
   Float_t         boosted250_lq0_dr;
   Float_t         boosted250_lq1_dr;
   Float_t         boosted250_qq_dr;
   Int_t           boosted500_decay_p4_;
   Float_t         boosted500_decay_p4_fCoordinates_fX[kMaxboosted500_decay_p4];   //[boosted500_decay_p4_]
   Float_t         boosted500_decay_p4_fCoordinates_fY[kMaxboosted500_decay_p4];   //[boosted500_decay_p4_]
   Float_t         boosted500_decay_p4_fCoordinates_fZ[kMaxboosted500_decay_p4];   //[boosted500_decay_p4_]
   Float_t         boosted500_decay_p4_fCoordinates_fT[kMaxboosted500_decay_p4];   //[boosted500_decay_p4_]
   vector<int>     *boosted500_decay_id;
   vector<int>     *boosted500_decay_isstar;
   vector<float>   *boosted500_decay_h_dr;
   vector<float>   *boosted500_decay_h_deta;
   vector<float>   *boosted500_decay_h_dphi;
   vector<float>   *boosted500_decay_h_deta_rotated;
   vector<float>   *boosted500_decay_h_dphi_rotated;
   Int_t           boosted500_lepton_p4_;
   Float_t         boosted500_lepton_p4_fCoordinates_fX[kMaxboosted500_lepton_p4];   //[boosted500_lepton_p4_]
   Float_t         boosted500_lepton_p4_fCoordinates_fY[kMaxboosted500_lepton_p4];   //[boosted500_lepton_p4_]
   Float_t         boosted500_lepton_p4_fCoordinates_fZ[kMaxboosted500_lepton_p4];   //[boosted500_lepton_p4_]
   Float_t         boosted500_lepton_p4_fCoordinates_fT[kMaxboosted500_lepton_p4];   //[boosted500_lepton_p4_]
   vector<int>     *boosted500_lepton_id;
   vector<int>     *boosted500_lepton_isstar;
   vector<float>   *boosted500_lepton_h_dr;
   vector<float>   *boosted500_lepton_h_deta;
   vector<float>   *boosted500_lepton_h_dphi;
   vector<float>   *boosted500_lepton_h_deta_rotated;
   vector<float>   *boosted500_lepton_h_dphi_rotated;
   Int_t           boosted500_quark_p4_;
   Float_t         boosted500_quark_p4_fCoordinates_fX[kMaxboosted500_quark_p4];   //[boosted500_quark_p4_]
   Float_t         boosted500_quark_p4_fCoordinates_fY[kMaxboosted500_quark_p4];   //[boosted500_quark_p4_]
   Float_t         boosted500_quark_p4_fCoordinates_fZ[kMaxboosted500_quark_p4];   //[boosted500_quark_p4_]
   Float_t         boosted500_quark_p4_fCoordinates_fT[kMaxboosted500_quark_p4];   //[boosted500_quark_p4_]
   vector<int>     *boosted500_quark_id;
   vector<int>     *boosted500_quark_isstar;
   vector<float>   *boosted500_quark_h_dr;
   vector<float>   *boosted500_quark_h_deta;
   vector<float>   *boosted500_quark_h_dphi;
   vector<float>   *boosted500_quark_h_deta_rotated;
   vector<float>   *boosted500_quark_h_dphi_rotated;
   Float_t         boosted500_lqq_max_dr;
   Float_t         boosted500_lq0_dr;
   Float_t         boosted500_lq1_dr;
   Float_t         boosted500_qq_dr;
   Int_t           boosted1000_decay_p4_;
   Float_t         boosted1000_decay_p4_fCoordinates_fX[kMaxboosted1000_decay_p4];   //[boosted1000_decay_p4_]
   Float_t         boosted1000_decay_p4_fCoordinates_fY[kMaxboosted1000_decay_p4];   //[boosted1000_decay_p4_]
   Float_t         boosted1000_decay_p4_fCoordinates_fZ[kMaxboosted1000_decay_p4];   //[boosted1000_decay_p4_]
   Float_t         boosted1000_decay_p4_fCoordinates_fT[kMaxboosted1000_decay_p4];   //[boosted1000_decay_p4_]
   vector<int>     *boosted1000_decay_id;
   vector<int>     *boosted1000_decay_isstar;
   vector<float>   *boosted1000_decay_h_dr;
   vector<float>   *boosted1000_decay_h_deta;
   vector<float>   *boosted1000_decay_h_dphi;
   vector<float>   *boosted1000_decay_h_deta_rotated;
   vector<float>   *boosted1000_decay_h_dphi_rotated;
   Int_t           boosted1000_lepton_p4_;
   Float_t         boosted1000_lepton_p4_fCoordinates_fX[kMaxboosted1000_lepton_p4];   //[boosted1000_lepton_p4_]
   Float_t         boosted1000_lepton_p4_fCoordinates_fY[kMaxboosted1000_lepton_p4];   //[boosted1000_lepton_p4_]
   Float_t         boosted1000_lepton_p4_fCoordinates_fZ[kMaxboosted1000_lepton_p4];   //[boosted1000_lepton_p4_]
   Float_t         boosted1000_lepton_p4_fCoordinates_fT[kMaxboosted1000_lepton_p4];   //[boosted1000_lepton_p4_]
   vector<int>     *boosted1000_lepton_id;
   vector<int>     *boosted1000_lepton_isstar;
   vector<float>   *boosted1000_lepton_h_dr;
   vector<float>   *boosted1000_lepton_h_deta;
   vector<float>   *boosted1000_lepton_h_dphi;
   vector<float>   *boosted1000_lepton_h_deta_rotated;
   vector<float>   *boosted1000_lepton_h_dphi_rotated;
   Int_t           boosted1000_quark_p4_;
   Float_t         boosted1000_quark_p4_fCoordinates_fX[kMaxboosted1000_quark_p4];   //[boosted1000_quark_p4_]
   Float_t         boosted1000_quark_p4_fCoordinates_fY[kMaxboosted1000_quark_p4];   //[boosted1000_quark_p4_]
   Float_t         boosted1000_quark_p4_fCoordinates_fZ[kMaxboosted1000_quark_p4];   //[boosted1000_quark_p4_]
   Float_t         boosted1000_quark_p4_fCoordinates_fT[kMaxboosted1000_quark_p4];   //[boosted1000_quark_p4_]
   vector<int>     *boosted1000_quark_id;
   vector<int>     *boosted1000_quark_isstar;
   vector<float>   *boosted1000_quark_h_dr;
   vector<float>   *boosted1000_quark_h_deta;
   vector<float>   *boosted1000_quark_h_dphi;
   vector<float>   *boosted1000_quark_h_deta_rotated;
   vector<float>   *boosted1000_quark_h_dphi_rotated;
   Float_t         boosted1000_lqq_max_dr;
   Float_t         boosted1000_lq0_dr;
   Float_t         boosted1000_lq1_dr;
   Float_t         boosted1000_qq_dr;
   Int_t           boosted1500_decay_p4_;
   Float_t         boosted1500_decay_p4_fCoordinates_fX[kMaxboosted1500_decay_p4];   //[boosted1500_decay_p4_]
   Float_t         boosted1500_decay_p4_fCoordinates_fY[kMaxboosted1500_decay_p4];   //[boosted1500_decay_p4_]
   Float_t         boosted1500_decay_p4_fCoordinates_fZ[kMaxboosted1500_decay_p4];   //[boosted1500_decay_p4_]
   Float_t         boosted1500_decay_p4_fCoordinates_fT[kMaxboosted1500_decay_p4];   //[boosted1500_decay_p4_]
   vector<int>     *boosted1500_decay_id;
   vector<int>     *boosted1500_decay_isstar;
   vector<float>   *boosted1500_decay_h_dr;
   vector<float>   *boosted1500_decay_h_deta;
   vector<float>   *boosted1500_decay_h_dphi;
   vector<float>   *boosted1500_decay_h_deta_rotated;
   vector<float>   *boosted1500_decay_h_dphi_rotated;
   Int_t           boosted1500_lepton_p4_;
   Float_t         boosted1500_lepton_p4_fCoordinates_fX[kMaxboosted1500_lepton_p4];   //[boosted1500_lepton_p4_]
   Float_t         boosted1500_lepton_p4_fCoordinates_fY[kMaxboosted1500_lepton_p4];   //[boosted1500_lepton_p4_]
   Float_t         boosted1500_lepton_p4_fCoordinates_fZ[kMaxboosted1500_lepton_p4];   //[boosted1500_lepton_p4_]
   Float_t         boosted1500_lepton_p4_fCoordinates_fT[kMaxboosted1500_lepton_p4];   //[boosted1500_lepton_p4_]
   vector<int>     *boosted1500_lepton_id;
   vector<int>     *boosted1500_lepton_isstar;
   vector<float>   *boosted1500_lepton_h_dr;
   vector<float>   *boosted1500_lepton_h_deta;
   vector<float>   *boosted1500_lepton_h_dphi;
   vector<float>   *boosted1500_lepton_h_deta_rotated;
   vector<float>   *boosted1500_lepton_h_dphi_rotated;
   Int_t           boosted1500_quark_p4_;
   Float_t         boosted1500_quark_p4_fCoordinates_fX[kMaxboosted1500_quark_p4];   //[boosted1500_quark_p4_]
   Float_t         boosted1500_quark_p4_fCoordinates_fY[kMaxboosted1500_quark_p4];   //[boosted1500_quark_p4_]
   Float_t         boosted1500_quark_p4_fCoordinates_fZ[kMaxboosted1500_quark_p4];   //[boosted1500_quark_p4_]
   Float_t         boosted1500_quark_p4_fCoordinates_fT[kMaxboosted1500_quark_p4];   //[boosted1500_quark_p4_]
   vector<int>     *boosted1500_quark_id;
   vector<int>     *boosted1500_quark_isstar;
   vector<float>   *boosted1500_quark_h_dr;
   vector<float>   *boosted1500_quark_h_deta;
   vector<float>   *boosted1500_quark_h_dphi;
   vector<float>   *boosted1500_quark_h_deta_rotated;
   vector<float>   *boosted1500_quark_h_dphi_rotated;
   Float_t         boosted1500_lqq_max_dr;
   Float_t         boosted1500_lq0_dr;
   Float_t         boosted1500_lq1_dr;
   Float_t         boosted1500_qq_dr;
   Int_t           iswhwww;
   Int_t           www_channel;
   Int_t           has_tau;
   Int_t           w_p4_;
   Float_t         w_p4_fCoordinates_fX[kMaxw_p4];   //[w_p4_]
   Float_t         w_p4_fCoordinates_fY[kMaxw_p4];   //[w_p4_]
   Float_t         w_p4_fCoordinates_fZ[kMaxw_p4];   //[w_p4_]
   Float_t         w_p4_fCoordinates_fT[kMaxw_p4];   //[w_p4_]
   vector<int>     *w_islep;
   vector<int>     *w_isstar;
   vector<int>     *w_isH;
   Int_t           l_p4_;
   Float_t         l_p4_fCoordinates_fX[kMaxl_p4];   //[l_p4_]
   Float_t         l_p4_fCoordinates_fY[kMaxl_p4];   //[l_p4_]
   Float_t         l_p4_fCoordinates_fZ[kMaxl_p4];   //[l_p4_]
   Float_t         l_p4_fCoordinates_fT[kMaxl_p4];   //[l_p4_]
   vector<float>   *l_w_pt;
   vector<float>   *l_w_eta;
   vector<float>   *l_w_phi;
   vector<float>   *l_w_mass;
   vector<int>     *l_w_id;
   vector<int>     *l_isstar;
   vector<int>     *l_isH;
   vector<int>     *l_istau;
   Int_t           q_p4_;
   Float_t         q_p4_fCoordinates_fX[kMaxq_p4];   //[q_p4_]
   Float_t         q_p4_fCoordinates_fY[kMaxq_p4];   //[q_p4_]
   Float_t         q_p4_fCoordinates_fZ[kMaxq_p4];   //[q_p4_]
   Float_t         q_p4_fCoordinates_fT[kMaxq_p4];   //[q_p4_]
   vector<float>   *q_w_pt;
   vector<float>   *q_w_eta;
   vector<float>   *q_w_phi;
   vector<float>   *q_w_mass;
   vector<int>     *q_w_id;
   vector<int>     *q_isstar;
   vector<int>     *q_isH;
   Float_t         dRllSS;
   Float_t         dRqqSS;
   Float_t         DPhill_higgs;
   Float_t         Mll_higgs;
   Float_t         MT_higgs;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_evt_scale1fb;   //!
   TBranch        *b_xsec_br;   //!
   TBranch        *b_evt_passgoodrunlist;   //!
   TBranch        *b_CMS4path;   //!
   TBranch        *b_CMS4index;   //!
   TBranch        *b_weight_fr_r1_f1;   //!
   TBranch        *b_weight_fr_r1_f2;   //!
   TBranch        *b_weight_fr_r1_f0p5;   //!
   TBranch        *b_weight_fr_r2_f1;   //!
   TBranch        *b_weight_fr_r2_f2;   //!
   TBranch        *b_weight_fr_r2_f0p5;   //!
   TBranch        *b_weight_fr_r0p5_f1;   //!
   TBranch        *b_weight_fr_r0p5_f2;   //!
   TBranch        *b_weight_fr_r0p5_f0p5;   //!
   TBranch        *b_weight_pdf_up;   //!
   TBranch        *b_weight_pdf_down;   //!
   TBranch        *b_weight_alphas_down;   //!
   TBranch        *b_weight_alphas_up;   //!
   TBranch        *b_weight_isr;   //!
   TBranch        *b_weight_isr_up;   //!
   TBranch        *b_weight_isr_down;   //!
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_DoubleEl_DZ;   //!
   TBranch        *b_HLT_DoubleEl_DZ_2;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_HLT_SingleEl8;   //!
   TBranch        *b_HLT_SingleEl17;   //!
   TBranch        *b_HLT_SingleIsoEl8;   //!
   TBranch        *b_HLT_SingleIsoEl17;   //!
   TBranch        *b_HLT_SingleIsoEl23;   //!
   TBranch        *b_HLT_SingleIsoMu8;   //!
   TBranch        *b_HLT_SingleIsoMu17;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_mc_HLT_DoubleMu;   //!
   TBranch        *b_mc_HLT_DoubleEl;   //!
   TBranch        *b_mc_HLT_DoubleEl_DZ;   //!
   TBranch        *b_mc_HLT_DoubleEl_DZ_2;   //!
   TBranch        *b_mc_HLT_MuEG;   //!
   TBranch        *b_mc_HLT_SingleEl8;   //!
   TBranch        *b_mc_HLT_SingleEl17;   //!
   TBranch        *b_mc_HLT_SingleIsoEl8;   //!
   TBranch        *b_mc_HLT_SingleIsoEl17;   //!
   TBranch        *b_mc_HLT_SingleIsoEl23;   //!
   TBranch        *b_mc_HLT_SingleIsoMu8;   //!
   TBranch        *b_mc_HLT_SingleIsoMu17;   //!
   TBranch        *b_mc_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_pass_duplicate_ee_em_mm;   //!
   TBranch        *b_pass_duplicate_mm_em_ee;   //!
   TBranch        *b_is2016;   //!
   TBranch        *b_is2017;   //!
   TBranch        *b_HLT_MuEG_2016;   //!
   TBranch        *b_mc_HLT_MuEG_2016;   //!
   TBranch        *b_pass_duplicate_ee_em2016_mm;   //!
   TBranch        *b_pass_duplicate_mm_em2016_ee;   //!
   TBranch        *b_passTrigger;   //!
   TBranch        *b_lep_p4_;   //!
   TBranch        *b_lep_p4_fCoordinates_fX;   //!
   TBranch        *b_lep_p4_fCoordinates_fY;   //!
   TBranch        *b_lep_p4_fCoordinates_fZ;   //!
   TBranch        *b_lep_p4_fCoordinates_fT;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_coneCorrPt;   //!
   TBranch        *b_lep_ip3d;   //!
   TBranch        *b_lep_ip3derr;   //!
   TBranch        *b_lep_isTriggerSafe_v1;   //!
   TBranch        *b_lep_lostHits;   //!
   TBranch        *b_lep_convVeto;   //!
   TBranch        *b_lep_motherIdSS;   //!
   TBranch        *b_lep_pass_VVV_cutbased_3l_fo;   //!
   TBranch        *b_lep_pass_VVV_cutbased_3l_tight;   //!
   TBranch        *b_lep_pass_VVV_cutbased_fo;   //!
   TBranch        *b_lep_pass_VVV_cutbased_tight;   //!
   TBranch        *b_lep_pass_VVV_cutbased_veto;   //!
   TBranch        *b_lep_pass_VVV_cutbased_fo_noiso;   //!
   TBranch        *b_lep_pass_VVV_cutbased_tight_noiso;   //!
   TBranch        *b_lep_pass_VVV_cutbased_veto_noiso;   //!
   TBranch        *b_lep_pass_POG_veto;   //!
   TBranch        *b_lep_pass_POG_loose;   //!
   TBranch        *b_lep_pass_POG_medium;   //!
   TBranch        *b_lep_pass_POG_tight;   //!
   TBranch        *b_lep_pdgId;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_pterr;   //!
   TBranch        *b_lep_relIso04DB;   //!
   TBranch        *b_lep_relIso03EA;   //!
   TBranch        *b_lep_relIso03EALep;   //!
   TBranch        *b_lep_relIso03EAv2;   //!
   TBranch        *b_lep_relIso04EAv2;   //!
   TBranch        *b_lep_relIso03EAv2Lep;   //!
   TBranch        *b_lep_tightCharge;   //!
   TBranch        *b_lep_trk_pt;   //!
   TBranch        *b_lep_charge;   //!
   TBranch        *b_lep_etaSC;   //!
   TBranch        *b_lep_MVA;   //!
   TBranch        *b_lep_isMediumPOG;   //!
   TBranch        *b_lep_isTightPOG;   //!
   TBranch        *b_lep_isFromW;   //!
   TBranch        *b_lep_isFromZ;   //!
   TBranch        *b_lep_isFromB;   //!
   TBranch        *b_lep_isFromC;   //!
   TBranch        *b_lep_isFromL;   //!
   TBranch        *b_lep_isFromLF;   //!
   TBranch        *b_lep_genPart_index;   //!
   TBranch        *b_lep_r9;   //!
   TBranch        *b_lep_nlayers;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_relIso03EA;   //!
   TBranch        *b_el_relIso03EALep;   //!
   TBranch        *b_el_ip3d;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_relIso04DB;   //!
   TBranch        *b_mu_relIso03EA;   //!
   TBranch        *b_mu_relIso03EALep;   //!
   TBranch        *b_mu_ip3d;   //!
   TBranch        *b_lbnt_pt;   //!
   TBranch        *b_lbnt_coneCorrPt;   //!
   TBranch        *b_lbnt_abseta;   //!
   TBranch        *b_lbnt_pdgId;   //!
   TBranch        *b_lbnt_el_pt;   //!
   TBranch        *b_lbnt_el_coneCorrPt;   //!
   TBranch        *b_lbnt_el_abseta;   //!
   TBranch        *b_lbnt_mu_pt;   //!
   TBranch        *b_lbnt_mu_coneCorrPt;   //!
   TBranch        *b_lbnt_mu_abseta;   //!
   TBranch        *b_jets_p4_;   //!
   TBranch        *b_jets_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_up_p4_;   //!
   TBranch        *b_jets_up_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_up_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_up_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_up_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_dn_p4_;   //!
   TBranch        *b_jets_dn_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_dn_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_dn_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_dn_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_csv;   //!
   TBranch        *b_jets_up_csv;   //!
   TBranch        *b_jets_dn_csv;   //!
   TBranch        *b_jets_jer_csv;   //!
   TBranch        *b_jets_jerup_csv;   //!
   TBranch        *b_jets_jerdn_csv;   //!
   TBranch        *b_jets_jer_p4_;   //!
   TBranch        *b_jets_jer_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_jer_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_jer_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_jer_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_jerup_p4_;   //!
   TBranch        *b_jets_jerup_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_jerup_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_jerup_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_jerup_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_jerdn_p4_;   //!
   TBranch        *b_jets_jerdn_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_jerdn_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_jerdn_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_jerdn_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_p4_;   //!
   TBranch        *b_jets30_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_up_p4_;   //!
   TBranch        *b_jets30_up_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_up_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_up_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_up_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_dn_p4_;   //!
   TBranch        *b_jets30_dn_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_dn_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_dn_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_dn_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_jer_p4_;   //!
   TBranch        *b_jets30_jer_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_jer_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_jer_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_jer_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_jerup_p4_;   //!
   TBranch        *b_jets30_jerup_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_jerup_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_jerup_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_jerup_p4_fCoordinates_fT;   //!
   TBranch        *b_jets30_jerdn_p4_;   //!
   TBranch        *b_jets30_jerdn_p4_fCoordinates_fX;   //!
   TBranch        *b_jets30_jerdn_p4_fCoordinates_fY;   //!
   TBranch        *b_jets30_jerdn_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets30_jerdn_p4_fCoordinates_fT;   //!
   TBranch        *b_ak8jets_p4_;   //!
   TBranch        *b_ak8jets_p4_fCoordinates_fX;   //!
   TBranch        *b_ak8jets_p4_fCoordinates_fY;   //!
   TBranch        *b_ak8jets_p4_fCoordinates_fZ;   //!
   TBranch        *b_ak8jets_p4_fCoordinates_fT;   //!
   TBranch        *b_ak8jets_softdropMass;   //!
   TBranch        *b_ak8jets_prunedMass;   //!
   TBranch        *b_ak8jets_trimmedMass;   //!
   TBranch        *b_ak8jets_mass;   //!
   TBranch        *b_ak8jets_nJettinessTau1;   //!
   TBranch        *b_ak8jets_nJettinessTau2;   //!
   TBranch        *b_ak8jets_softdropPuppiSubjet1;   //!
   TBranch        *b_ak8jets_softdropPuppiSubjet2;   //!
   TBranch        *b_ak8jets_puppi_softdropMass;   //!
   TBranch        *b_ak8jets_puppi_nJettinessTau1;   //!
   TBranch        *b_ak8jets_puppi_nJettinessTau2;   //!
   TBranch        *b_ak8jets_puppi_eta;   //!
   TBranch        *b_ak8jets_puppi_phi;   //!
   TBranch        *b_ak8jets_puppi_pt;   //!
   TBranch        *b_ak8jets_puppi_mass;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_up_pt;   //!
   TBranch        *b_met_up_phi;   //!
   TBranch        *b_met_dn_pt;   //!
   TBranch        *b_met_dn_phi;   //!
   TBranch        *b_met_gen_pt;   //!
   TBranch        *b_met_gen_phi;   //!
   TBranch        *b_met_jer_pt;   //!
   TBranch        *b_met_jerup_pt;   //!
   TBranch        *b_met_jerdn_pt;   //!
   TBranch        *b_met_jer_phi;   //!
   TBranch        *b_met_jerup_phi;   //!
   TBranch        *b_met_jerdn_phi;   //!
   TBranch        *b_firstgoodvertex;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_nVert;   //!
   TBranch        *b_nisoTrack_mt2_cleaned_VVV_cutbased_veto;   //!
   TBranch        *b_weight_btagsf;   //!
   TBranch        *b_weight_btagsf_heavy_DN;   //!
   TBranch        *b_weight_btagsf_heavy_UP;   //!
   TBranch        *b_weight_btagsf_light_DN;   //!
   TBranch        *b_weight_btagsf_light_UP;   //!
   TBranch        *b_gen_ht;   //!
   TBranch        *b_genPart_p4_;   //!
   TBranch        *b_genPart_p4_fCoordinates_fX;   //!
   TBranch        *b_genPart_p4_fCoordinates_fY;   //!
   TBranch        *b_genPart_p4_fCoordinates_fZ;   //!
   TBranch        *b_genPart_p4_fCoordinates_fT;   //!
   TBranch        *b_genPart_motherId;   //!
   TBranch        *b_genPart_pdgId;   //!
   TBranch        *b_genPart_charge;   //!
   TBranch        *b_genPart_status;   //!
   TBranch        *b_ngenLep;   //!
   TBranch        *b_ngenLepFromTau;   //!
   TBranch        *b_ngenLepFromBoson;   //!
   TBranch        *b_Flag_AllEventFilters;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_HBHEIsoNoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_badChargedCandidateFilter;   //!
   TBranch        *b_Flag_badMuonFilter;   //!
   TBranch        *b_Flag_badMuonFilterv2;   //!
   TBranch        *b_Flag_badChargedCandidateFilterv2;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_globalTightHalo2016;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_trackingFailureFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_badMuons;   //!
   TBranch        *b_Flag_duplicateMuons;   //!
   TBranch        *b_Flag_noBadMuons;   //!
   TBranch        *b_fastsimfilt;   //!
   TBranch        *b_nVlep;   //!
   TBranch        *b_nTlep;   //!
   TBranch        *b_nTlepSS;   //!
   TBranch        *b_nLlep;   //!
   TBranch        *b_nLlep3L;   //!
   TBranch        *b_nTlep3L;   //!
   TBranch        *b_nSFOS;   //!
   TBranch        *b_nSFOSinZ;   //!
   TBranch        *b_nj;   //!
   TBranch        *b_nj_up;   //!
   TBranch        *b_nj_dn;   //!
   TBranch        *b_nj_jer;   //!
   TBranch        *b_nj_jerup;   //!
   TBranch        *b_nj_jerdn;   //!
   TBranch        *b_nj30;   //!
   TBranch        *b_nj30_up;   //!
   TBranch        *b_nj30_dn;   //!
   TBranch        *b_nj30_jer;   //!
   TBranch        *b_nj30_jerup;   //!
   TBranch        *b_nj30_jerdn;   //!
   TBranch        *b_nb;   //!
   TBranch        *b_nb_up;   //!
   TBranch        *b_nb_dn;   //!
   TBranch        *b_nb_jer;   //!
   TBranch        *b_nb_jerup;   //!
   TBranch        *b_nb_jerdn;   //!
   TBranch        *b_Ml0j0;   //!
   TBranch        *b_Ml0j0_up;   //!
   TBranch        *b_Ml0j0_dn;   //!
   TBranch        *b_Ml0j0_jer;   //!
   TBranch        *b_Ml0j0_jerup;   //!
   TBranch        *b_Ml0j0_jerdn;   //!
   TBranch        *b_Ml0j1;   //!
   TBranch        *b_Ml0j1_up;   //!
   TBranch        *b_Ml0j1_dn;   //!
   TBranch        *b_Ml0j1_jer;   //!
   TBranch        *b_Ml0j1_jerup;   //!
   TBranch        *b_Ml0j1_jerdn;   //!
   TBranch        *b_Ml1j0;   //!
   TBranch        *b_Ml1j0_up;   //!
   TBranch        *b_Ml1j0_dn;   //!
   TBranch        *b_Ml1j0_jer;   //!
   TBranch        *b_Ml1j0_jerup;   //!
   TBranch        *b_Ml1j0_jerdn;   //!
   TBranch        *b_Ml1j1;   //!
   TBranch        *b_Ml1j1_up;   //!
   TBranch        *b_Ml1j1_dn;   //!
   TBranch        *b_Ml1j1_jer;   //!
   TBranch        *b_Ml1j1_jerup;   //!
   TBranch        *b_Ml1j1_jerdn;   //!
   TBranch        *b_MinMlj;   //!
   TBranch        *b_MinMlj_up;   //!
   TBranch        *b_MinMlj_dn;   //!
   TBranch        *b_MinMlj_jer;   //!
   TBranch        *b_MinMlj_jerup;   //!
   TBranch        *b_MinMlj_jerdn;   //!
   TBranch        *b_SumMinMlj01;   //!
   TBranch        *b_SumMinMlj01_up;   //!
   TBranch        *b_SumMinMlj01_dn;   //!
   TBranch        *b_SumMinMlj01_jer;   //!
   TBranch        *b_SumMinMlj01_jerup;   //!
   TBranch        *b_SumMinMlj01_jerdn;   //!
   TBranch        *b_MaxMlj;   //!
   TBranch        *b_MaxMlj_up;   //!
   TBranch        *b_MaxMlj_dn;   //!
   TBranch        *b_MaxMlj_jer;   //!
   TBranch        *b_MaxMlj_jerup;   //!
   TBranch        *b_MaxMlj_jerdn;   //!
   TBranch        *b_SumMlj;   //!
   TBranch        *b_SumMlj_up;   //!
   TBranch        *b_SumMlj_dn;   //!
   TBranch        *b_SumMlj_jer;   //!
   TBranch        *b_SumMlj_jerup;   //!
   TBranch        *b_SumMlj_jerdn;   //!
   TBranch        *b_Ml0jj;   //!
   TBranch        *b_Ml0jj_up;   //!
   TBranch        *b_Ml0jj_dn;   //!
   TBranch        *b_Ml0jj_jer;   //!
   TBranch        *b_Ml0jj_jerup;   //!
   TBranch        *b_Ml0jj_jerdn;   //!
   TBranch        *b_Ml1jj;   //!
   TBranch        *b_Ml1jj_up;   //!
   TBranch        *b_Ml1jj_dn;   //!
   TBranch        *b_Ml1jj_jer;   //!
   TBranch        *b_Ml1jj_jerup;   //!
   TBranch        *b_Ml1jj_jerdn;   //!
   TBranch        *b_MinMljj;   //!
   TBranch        *b_MinMljj_up;   //!
   TBranch        *b_MinMljj_dn;   //!
   TBranch        *b_MinMljj_jer;   //!
   TBranch        *b_MinMljj_jerup;   //!
   TBranch        *b_MinMljj_jerdn;   //!
   TBranch        *b_MaxMljj;   //!
   TBranch        *b_MaxMljj_up;   //!
   TBranch        *b_MaxMljj_dn;   //!
   TBranch        *b_MaxMljj_jer;   //!
   TBranch        *b_MaxMljj_jerup;   //!
   TBranch        *b_MaxMljj_jerdn;   //!
   TBranch        *b_SumMljj;   //!
   TBranch        *b_SumMljj_up;   //!
   TBranch        *b_SumMljj_dn;   //!
   TBranch        *b_SumMljj_jer;   //!
   TBranch        *b_SumMljj_jerup;   //!
   TBranch        *b_SumMljj_jerdn;   //!
   TBranch        *b_Mjj;   //!
   TBranch        *b_Mjj_up;   //!
   TBranch        *b_Mjj_dn;   //!
   TBranch        *b_Mjj_jer;   //!
   TBranch        *b_Mjj_jerup;   //!
   TBranch        *b_Mjj_jerdn;   //!
   TBranch        *b_DRjj;   //!
   TBranch        *b_DRjj_up;   //!
   TBranch        *b_DRjj_dn;   //!
   TBranch        *b_DRjj_jer;   //!
   TBranch        *b_DRjj_jerup;   //!
   TBranch        *b_DRjj_jerdn;   //!
   TBranch        *b_jet0_wtag_p4_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_up_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_up_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_up_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_up_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_dn_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_dn_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_dn_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_dn_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_jer_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_jer_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_jer_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_jer_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_jerup_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_jerup_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_jerup_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_jerup_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_jerdn_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_jerdn_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_jerdn_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_jerdn_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_up_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_up_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_up_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_up_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_dn_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_dn_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_dn_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_dn_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_jer_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_jer_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_jer_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_jer_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_jerup_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_jerup_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_jerup_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_jerup_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_jerdn_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_jerdn_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_jerdn_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_jerdn_fCoordinates_fT;   //!
   TBranch        *b_MjjDR1;   //!
   TBranch        *b_MjjDR1_up;   //!
   TBranch        *b_MjjDR1_dn;   //!
   TBranch        *b_MjjDR1_jer;   //!
   TBranch        *b_MjjDR1_jerup;   //!
   TBranch        *b_MjjDR1_jerdn;   //!
   TBranch        *b_DRjjDR1;   //!
   TBranch        *b_DRjjDR1_up;   //!
   TBranch        *b_DRjjDR1_dn;   //!
   TBranch        *b_DRjjDR1_jer;   //!
   TBranch        *b_DRjjDR1_jerup;   //!
   TBranch        *b_DRjjDR1_jerdn;   //!
   TBranch        *b_jet0_wtag_p4_DR1_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_DR1_up_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_up_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_up_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_up_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_DR1_dn_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_dn_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_dn_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_dn_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jer_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jer_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jer_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jer_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerup_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerup_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerup_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerup_fCoordinates_fT;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fX;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fY;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fZ;   //!
   TBranch        *b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_up_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_up_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_up_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_up_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_dn_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_dn_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_dn_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_dn_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jer_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jer_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jer_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jer_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerup_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerup_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerup_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerup_fCoordinates_fT;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fX;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fY;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fZ;   //!
   TBranch        *b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fT;   //!
   TBranch        *b_MjjVBF;   //!
   TBranch        *b_MjjVBF_up;   //!
   TBranch        *b_MjjVBF_dn;   //!
   TBranch        *b_MjjVBF_jer;   //!
   TBranch        *b_MjjVBF_jerup;   //!
   TBranch        *b_MjjVBF_jerdn;   //!
   TBranch        *b_DetajjVBF;   //!
   TBranch        *b_DetajjVBF_up;   //!
   TBranch        *b_DetajjVBF_dn;   //!
   TBranch        *b_DetajjVBF_jer;   //!
   TBranch        *b_DetajjVBF_jerup;   //!
   TBranch        *b_DetajjVBF_jerdn;   //!
   TBranch        *b_MjjL;   //!
   TBranch        *b_MjjL_up;   //!
   TBranch        *b_MjjL_dn;   //!
   TBranch        *b_MjjL_jer;   //!
   TBranch        *b_MjjL_jerup;   //!
   TBranch        *b_MjjL_jerdn;   //!
   TBranch        *b_DetajjL;   //!
   TBranch        *b_DetajjL_up;   //!
   TBranch        *b_DetajjL_dn;   //!
   TBranch        *b_DetajjL_jer;   //!
   TBranch        *b_DetajjL_jerup;   //!
   TBranch        *b_DetajjL_jerdn;   //!
   TBranch        *b_MllSS;   //!
   TBranch        *b_MeeSS;   //!
   TBranch        *b_Mll3L;   //!
   TBranch        *b_Mee3L;   //!
   TBranch        *b_Mll3L1;   //!
   TBranch        *b_M3l;   //!
   TBranch        *b_Pt3l;   //!
   TBranch        *b_M01;   //!
   TBranch        *b_M02;   //!
   TBranch        *b_M12;   //!
   TBranch        *b_isSFOS01;   //!
   TBranch        *b_isSFOS02;   //!
   TBranch        *b_isSFOS12;   //!
   TBranch        *b_DPhi3lMET;   //!
   TBranch        *b_DPhi3lMET_up;   //!
   TBranch        *b_DPhi3lMET_dn;   //!
   TBranch        *b_DPhi3lMET_jer;   //!
   TBranch        *b_DPhi3lMET_jerup;   //!
   TBranch        *b_DPhi3lMET_jerdn;   //!
   TBranch        *b_DPhi3lMET_gen;   //!
   TBranch        *b_MTmax;   //!
   TBranch        *b_MTmax_up;   //!
   TBranch        *b_MTmax_dn;   //!
   TBranch        *b_MTmax_jer;   //!
   TBranch        *b_MTmax_jerup;   //!
   TBranch        *b_MTmax_jerdn;   //!
   TBranch        *b_MTmax_gen;   //!
   TBranch        *b_MTmin;   //!
   TBranch        *b_MTmin_up;   //!
   TBranch        *b_MTmin_dn;   //!
   TBranch        *b_MTmin_jer;   //!
   TBranch        *b_MTmin_jerup;   //!
   TBranch        *b_MTmin_jerdn;   //!
   TBranch        *b_MTmin_gen;   //!
   TBranch        *b_MT3rd;   //!
   TBranch        *b_MT3rd_up;   //!
   TBranch        *b_MT3rd_dn;   //!
   TBranch        *b_MT3rd_jer;   //!
   TBranch        *b_MT3rd_jerup;   //!
   TBranch        *b_MT3rd_jerdn;   //!
   TBranch        *b_MT3rd_gen;   //!
   TBranch        *b_MTmax3L;   //!
   TBranch        *b_MTmax3L_up;   //!
   TBranch        *b_MTmax3L_dn;   //!
   TBranch        *b_MTmax3L_jer;   //!
   TBranch        *b_MTmax3L_jerup;   //!
   TBranch        *b_MTmax3L_jerdn;   //!
   TBranch        *b_MTmax3L_gen;   //!
   TBranch        *b_passSSee;   //!
   TBranch        *b_passSSem;   //!
   TBranch        *b_passSSmm;   //!
   TBranch        *b_lep_idx0_SS;   //!
   TBranch        *b_lep_idx1_SS;   //!
   TBranch        *b_bkgtype;   //!
   TBranch        *b_vetophoton;   //!
   TBranch        *b_purewgt;   //!
   TBranch        *b_purewgt_up;   //!
   TBranch        *b_purewgt_dn;   //!
   TBranch        *b_ffwgt;   //!
   TBranch        *b_ffwgt_up;   //!
   TBranch        *b_ffwgt_dn;   //!
   TBranch        *b_ffwgt_el_up;   //!
   TBranch        *b_ffwgt_el_dn;   //!
   TBranch        *b_ffwgt_mu_up;   //!
   TBranch        *b_ffwgt_mu_dn;   //!
   TBranch        *b_ffwgt_closure_up;   //!
   TBranch        *b_ffwgt_closure_dn;   //!
   TBranch        *b_ffwgt_closure_el_up;   //!
   TBranch        *b_ffwgt_closure_el_dn;   //!
   TBranch        *b_ffwgt_closure_mu_up;   //!
   TBranch        *b_ffwgt_closure_mu_dn;   //!
   TBranch        *b_ffwgt_full_up;   //!
   TBranch        *b_ffwgt_full_dn;   //!
   TBranch        *b_ffwgtqcd;   //!
   TBranch        *b_ffwgtqcd_up;   //!
   TBranch        *b_ffwgtqcd_dn;   //!
   TBranch        *b_lepsf;   //!
   TBranch        *b_lepsf_up;   //!
   TBranch        *b_lepsf_dn;   //!
   TBranch        *b_trigeff;   //!
   TBranch        *b_trigeff_up;   //!
   TBranch        *b_trigeff_dn;   //!
   TBranch        *b_trigsf;   //!
   TBranch        *b_trigsf_up;   //!
   TBranch        *b_trigsf_dn;   //!
   TBranch        *b_musmear_sf;   //!
   TBranch        *b_higgsdecay;   //!
   TBranch        *b_nlep;   //!
   TBranch        *b_nquark;   //!
   TBranch        *b_wa_id;   //!
   TBranch        *b_higgs_p4_fCoordinates_fX;   //!
   TBranch        *b_higgs_p4_fCoordinates_fY;   //!
   TBranch        *b_higgs_p4_fCoordinates_fZ;   //!
   TBranch        *b_higgs_p4_fCoordinates_fT;   //!
   TBranch        *b_decay_p4_;   //!
   TBranch        *b_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_decay_id;   //!
   TBranch        *b_decay_isstar;   //!
   TBranch        *b_lepton_p4_;   //!
   TBranch        *b_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_lepton_id;   //!
   TBranch        *b_lepton_isstar;   //!
   TBranch        *b_quark_p4_;   //!
   TBranch        *b_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_quark_id;   //!
   TBranch        *b_quark_isstar;   //!
   TBranch        *b_lqq_max_dr;   //!
   TBranch        *b_lq0_dr;   //!
   TBranch        *b_lq1_dr;   //!
   TBranch        *b_qq_dr;   //!
   TBranch        *b_boosted0_decay_p4_;   //!
   TBranch        *b_boosted0_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted0_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted0_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted0_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted0_decay_id;   //!
   TBranch        *b_boosted0_decay_isstar;   //!
   TBranch        *b_boosted0_decay_h_dr;   //!
   TBranch        *b_boosted0_decay_h_deta;   //!
   TBranch        *b_boosted0_decay_h_dphi;   //!
   TBranch        *b_boosted0_decay_h_deta_rotated;   //!
   TBranch        *b_boosted0_decay_h_dphi_rotated;   //!
   TBranch        *b_boosted0_lepton_p4_;   //!
   TBranch        *b_boosted0_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted0_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted0_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted0_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted0_lepton_id;   //!
   TBranch        *b_boosted0_lepton_isstar;   //!
   TBranch        *b_boosted0_lepton_h_dr;   //!
   TBranch        *b_boosted0_lepton_h_deta;   //!
   TBranch        *b_boosted0_lepton_h_dphi;   //!
   TBranch        *b_boosted0_lepton_h_deta_rotated;   //!
   TBranch        *b_boosted0_lepton_h_dphi_rotated;   //!
   TBranch        *b_boosted0_quark_p4_;   //!
   TBranch        *b_boosted0_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted0_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted0_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted0_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted0_quark_id;   //!
   TBranch        *b_boosted0_quark_isstar;   //!
   TBranch        *b_boosted0_quark_h_dr;   //!
   TBranch        *b_boosted0_quark_h_deta;   //!
   TBranch        *b_boosted0_quark_h_dphi;   //!
   TBranch        *b_boosted0_quark_h_deta_rotated;   //!
   TBranch        *b_boosted0_quark_h_dphi_rotated;   //!
   TBranch        *b_boosted250_decay_p4_;   //!
   TBranch        *b_boosted250_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted250_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted250_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted250_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted250_decay_id;   //!
   TBranch        *b_boosted250_decay_isstar;   //!
   TBranch        *b_boosted250_decay_h_dr;   //!
   TBranch        *b_boosted250_decay_h_deta;   //!
   TBranch        *b_boosted250_decay_h_dphi;   //!
   TBranch        *b_boosted250_decay_h_deta_rotated;   //!
   TBranch        *b_boosted250_decay_h_dphi_rotated;   //!
   TBranch        *b_boosted250_lepton_p4_;   //!
   TBranch        *b_boosted250_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted250_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted250_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted250_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted250_lepton_id;   //!
   TBranch        *b_boosted250_lepton_isstar;   //!
   TBranch        *b_boosted250_lepton_h_dr;   //!
   TBranch        *b_boosted250_lepton_h_deta;   //!
   TBranch        *b_boosted250_lepton_h_dphi;   //!
   TBranch        *b_boosted250_lepton_h_deta_rotated;   //!
   TBranch        *b_boosted250_lepton_h_dphi_rotated;   //!
   TBranch        *b_boosted250_quark_p4_;   //!
   TBranch        *b_boosted250_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted250_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted250_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted250_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted250_quark_id;   //!
   TBranch        *b_boosted250_quark_isstar;   //!
   TBranch        *b_boosted250_quark_h_dr;   //!
   TBranch        *b_boosted250_quark_h_deta;   //!
   TBranch        *b_boosted250_quark_h_dphi;   //!
   TBranch        *b_boosted250_quark_h_deta_rotated;   //!
   TBranch        *b_boosted250_quark_h_dphi_rotated;   //!
   TBranch        *b_boosted250_lqq_max_dr;   //!
   TBranch        *b_boosted250_lq0_dr;   //!
   TBranch        *b_boosted250_lq1_dr;   //!
   TBranch        *b_boosted250_qq_dr;   //!
   TBranch        *b_boosted500_decay_p4_;   //!
   TBranch        *b_boosted500_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted500_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted500_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted500_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted500_decay_id;   //!
   TBranch        *b_boosted500_decay_isstar;   //!
   TBranch        *b_boosted500_decay_h_dr;   //!
   TBranch        *b_boosted500_decay_h_deta;   //!
   TBranch        *b_boosted500_decay_h_dphi;   //!
   TBranch        *b_boosted500_decay_h_deta_rotated;   //!
   TBranch        *b_boosted500_decay_h_dphi_rotated;   //!
   TBranch        *b_boosted500_lepton_p4_;   //!
   TBranch        *b_boosted500_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted500_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted500_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted500_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted500_lepton_id;   //!
   TBranch        *b_boosted500_lepton_isstar;   //!
   TBranch        *b_boosted500_lepton_h_dr;   //!
   TBranch        *b_boosted500_lepton_h_deta;   //!
   TBranch        *b_boosted500_lepton_h_dphi;   //!
   TBranch        *b_boosted500_lepton_h_deta_rotated;   //!
   TBranch        *b_boosted500_lepton_h_dphi_rotated;   //!
   TBranch        *b_boosted500_quark_p4_;   //!
   TBranch        *b_boosted500_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted500_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted500_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted500_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted500_quark_id;   //!
   TBranch        *b_boosted500_quark_isstar;   //!
   TBranch        *b_boosted500_quark_h_dr;   //!
   TBranch        *b_boosted500_quark_h_deta;   //!
   TBranch        *b_boosted500_quark_h_dphi;   //!
   TBranch        *b_boosted500_quark_h_deta_rotated;   //!
   TBranch        *b_boosted500_quark_h_dphi_rotated;   //!
   TBranch        *b_boosted500_lqq_max_dr;   //!
   TBranch        *b_boosted500_lq0_dr;   //!
   TBranch        *b_boosted500_lq1_dr;   //!
   TBranch        *b_boosted500_qq_dr;   //!
   TBranch        *b_boosted1000_decay_p4_;   //!
   TBranch        *b_boosted1000_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1000_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1000_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1000_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1000_decay_id;   //!
   TBranch        *b_boosted1000_decay_isstar;   //!
   TBranch        *b_boosted1000_decay_h_dr;   //!
   TBranch        *b_boosted1000_decay_h_deta;   //!
   TBranch        *b_boosted1000_decay_h_dphi;   //!
   TBranch        *b_boosted1000_decay_h_deta_rotated;   //!
   TBranch        *b_boosted1000_decay_h_dphi_rotated;   //!
   TBranch        *b_boosted1000_lepton_p4_;   //!
   TBranch        *b_boosted1000_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1000_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1000_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1000_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1000_lepton_id;   //!
   TBranch        *b_boosted1000_lepton_isstar;   //!
   TBranch        *b_boosted1000_lepton_h_dr;   //!
   TBranch        *b_boosted1000_lepton_h_deta;   //!
   TBranch        *b_boosted1000_lepton_h_dphi;   //!
   TBranch        *b_boosted1000_lepton_h_deta_rotated;   //!
   TBranch        *b_boosted1000_lepton_h_dphi_rotated;   //!
   TBranch        *b_boosted1000_quark_p4_;   //!
   TBranch        *b_boosted1000_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1000_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1000_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1000_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1000_quark_id;   //!
   TBranch        *b_boosted1000_quark_isstar;   //!
   TBranch        *b_boosted1000_quark_h_dr;   //!
   TBranch        *b_boosted1000_quark_h_deta;   //!
   TBranch        *b_boosted1000_quark_h_dphi;   //!
   TBranch        *b_boosted1000_quark_h_deta_rotated;   //!
   TBranch        *b_boosted1000_quark_h_dphi_rotated;   //!
   TBranch        *b_boosted1000_lqq_max_dr;   //!
   TBranch        *b_boosted1000_lq0_dr;   //!
   TBranch        *b_boosted1000_lq1_dr;   //!
   TBranch        *b_boosted1000_qq_dr;   //!
   TBranch        *b_boosted1500_decay_p4_;   //!
   TBranch        *b_boosted1500_decay_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1500_decay_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1500_decay_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1500_decay_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1500_decay_id;   //!
   TBranch        *b_boosted1500_decay_isstar;   //!
   TBranch        *b_boosted1500_decay_h_dr;   //!
   TBranch        *b_boosted1500_decay_h_deta;   //!
   TBranch        *b_boosted1500_decay_h_dphi;   //!
   TBranch        *b_boosted1500_decay_h_deta_rotated;   //!
   TBranch        *b_boosted1500_decay_h_dphi_rotated;   //!
   TBranch        *b_boosted1500_lepton_p4_;   //!
   TBranch        *b_boosted1500_lepton_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1500_lepton_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1500_lepton_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1500_lepton_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1500_lepton_id;   //!
   TBranch        *b_boosted1500_lepton_isstar;   //!
   TBranch        *b_boosted1500_lepton_h_dr;   //!
   TBranch        *b_boosted1500_lepton_h_deta;   //!
   TBranch        *b_boosted1500_lepton_h_dphi;   //!
   TBranch        *b_boosted1500_lepton_h_deta_rotated;   //!
   TBranch        *b_boosted1500_lepton_h_dphi_rotated;   //!
   TBranch        *b_boosted1500_quark_p4_;   //!
   TBranch        *b_boosted1500_quark_p4_fCoordinates_fX;   //!
   TBranch        *b_boosted1500_quark_p4_fCoordinates_fY;   //!
   TBranch        *b_boosted1500_quark_p4_fCoordinates_fZ;   //!
   TBranch        *b_boosted1500_quark_p4_fCoordinates_fT;   //!
   TBranch        *b_boosted1500_quark_id;   //!
   TBranch        *b_boosted1500_quark_isstar;   //!
   TBranch        *b_boosted1500_quark_h_dr;   //!
   TBranch        *b_boosted1500_quark_h_deta;   //!
   TBranch        *b_boosted1500_quark_h_dphi;   //!
   TBranch        *b_boosted1500_quark_h_deta_rotated;   //!
   TBranch        *b_boosted1500_quark_h_dphi_rotated;   //!
   TBranch        *b_boosted1500_lqq_max_dr;   //!
   TBranch        *b_boosted1500_lq0_dr;   //!
   TBranch        *b_boosted1500_lq1_dr;   //!
   TBranch        *b_boosted1500_qq_dr;   //!
   TBranch        *b_iswhwww;   //!
   TBranch        *b_www_channel;   //!
   TBranch        *b_has_tau;   //!
   TBranch        *b_w_p4_;   //!
   TBranch        *b_w_p4_fCoordinates_fX;   //!
   TBranch        *b_w_p4_fCoordinates_fY;   //!
   TBranch        *b_w_p4_fCoordinates_fZ;   //!
   TBranch        *b_w_p4_fCoordinates_fT;   //!
   TBranch        *b_w_islep;   //!
   TBranch        *b_w_isstar;   //!
   TBranch        *b_w_isH;   //!
   TBranch        *b_l_p4_;   //!
   TBranch        *b_l_p4_fCoordinates_fX;   //!
   TBranch        *b_l_p4_fCoordinates_fY;   //!
   TBranch        *b_l_p4_fCoordinates_fZ;   //!
   TBranch        *b_l_p4_fCoordinates_fT;   //!
   TBranch        *b_l_w_pt;   //!
   TBranch        *b_l_w_eta;   //!
   TBranch        *b_l_w_phi;   //!
   TBranch        *b_l_w_mass;   //!
   TBranch        *b_l_w_id;   //!
   TBranch        *b_l_isstar;   //!
   TBranch        *b_l_isH;   //!
   TBranch        *b_l_istau;   //!
   TBranch        *b_q_p4_;   //!
   TBranch        *b_q_p4_fCoordinates_fX;   //!
   TBranch        *b_q_p4_fCoordinates_fY;   //!
   TBranch        *b_q_p4_fCoordinates_fZ;   //!
   TBranch        *b_q_p4_fCoordinates_fT;   //!
   TBranch        *b_q_w_pt;   //!
   TBranch        *b_q_w_eta;   //!
   TBranch        *b_q_w_phi;   //!
   TBranch        *b_q_w_mass;   //!
   TBranch        *b_q_w_id;   //!
   TBranch        *b_q_isstar;   //!
   TBranch        *b_q_isH;   //!
   TBranch        *b_dRllSS;   //!
   TBranch        *b_dRqqSS;   //!
   TBranch        *b_DPhill_higgs;   //!
   TBranch        *b_Mll_higgs;   //!
   TBranch        *b_MT_higgs;   //!

   t_www(TTree *tree=0);
   virtual ~t_www();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t_www_cxx
t_www::t_www(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("results.root");
      }
      f->GetObject("t_www",tree);

   }
   Init(tree);
}

t_www::~t_www()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t_www::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t_www::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void t_www::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   CMS4path = 0;
   lep_pt = 0;
   lep_eta = 0;
   lep_phi = 0;
   lep_coneCorrPt = 0;
   lep_ip3d = 0;
   lep_ip3derr = 0;
   lep_isTriggerSafe_v1 = 0;
   lep_lostHits = 0;
   lep_convVeto = 0;
   lep_motherIdSS = 0;
   lep_pass_VVV_cutbased_3l_fo = 0;
   lep_pass_VVV_cutbased_3l_tight = 0;
   lep_pass_VVV_cutbased_fo = 0;
   lep_pass_VVV_cutbased_tight = 0;
   lep_pass_VVV_cutbased_veto = 0;
   lep_pass_VVV_cutbased_fo_noiso = 0;
   lep_pass_VVV_cutbased_tight_noiso = 0;
   lep_pass_VVV_cutbased_veto_noiso = 0;
   lep_pass_POG_veto = 0;
   lep_pass_POG_loose = 0;
   lep_pass_POG_medium = 0;
   lep_pass_POG_tight = 0;
   lep_pdgId = 0;
   lep_dxy = 0;
   lep_dz = 0;
   lep_pterr = 0;
   lep_relIso04DB = 0;
   lep_relIso03EA = 0;
   lep_relIso03EALep = 0;
   lep_relIso03EAv2 = 0;
   lep_relIso04EAv2 = 0;
   lep_relIso03EAv2Lep = 0;
   lep_tightCharge = 0;
   lep_trk_pt = 0;
   lep_charge = 0;
   lep_etaSC = 0;
   lep_MVA = 0;
   lep_isMediumPOG = 0;
   lep_isTightPOG = 0;
   lep_isFromW = 0;
   lep_isFromZ = 0;
   lep_isFromB = 0;
   lep_isFromC = 0;
   lep_isFromL = 0;
   lep_isFromLF = 0;
   lep_genPart_index = 0;
   lep_r9 = 0;
   lep_nlayers = 0;
   jets_csv = 0;
   jets_up_csv = 0;
   jets_dn_csv = 0;
   jets_jer_csv = 0;
   jets_jerup_csv = 0;
   jets_jerdn_csv = 0;
   ak8jets_softdropMass = 0;
   ak8jets_prunedMass = 0;
   ak8jets_trimmedMass = 0;
   ak8jets_mass = 0;
   ak8jets_nJettinessTau1 = 0;
   ak8jets_nJettinessTau2 = 0;
   ak8jets_softdropPuppiSubjet1 = 0;
   ak8jets_softdropPuppiSubjet2 = 0;
   ak8jets_puppi_softdropMass = 0;
   ak8jets_puppi_nJettinessTau1 = 0;
   ak8jets_puppi_nJettinessTau2 = 0;
   ak8jets_puppi_eta = 0;
   ak8jets_puppi_phi = 0;
   ak8jets_puppi_pt = 0;
   ak8jets_puppi_mass = 0;
   genPart_motherId = 0;
   genPart_pdgId = 0;
   genPart_charge = 0;
   genPart_status = 0;
   bkgtype = 0;
   decay_id = 0;
   decay_isstar = 0;
   lepton_id = 0;
   lepton_isstar = 0;
   quark_id = 0;
   quark_isstar = 0;
   boosted0_decay_id = 0;
   boosted0_decay_isstar = 0;
   boosted0_decay_h_dr = 0;
   boosted0_decay_h_deta = 0;
   boosted0_decay_h_dphi = 0;
   boosted0_decay_h_deta_rotated = 0;
   boosted0_decay_h_dphi_rotated = 0;
   boosted0_lepton_id = 0;
   boosted0_lepton_isstar = 0;
   boosted0_lepton_h_dr = 0;
   boosted0_lepton_h_deta = 0;
   boosted0_lepton_h_dphi = 0;
   boosted0_lepton_h_deta_rotated = 0;
   boosted0_lepton_h_dphi_rotated = 0;
   boosted0_quark_id = 0;
   boosted0_quark_isstar = 0;
   boosted0_quark_h_dr = 0;
   boosted0_quark_h_deta = 0;
   boosted0_quark_h_dphi = 0;
   boosted0_quark_h_deta_rotated = 0;
   boosted0_quark_h_dphi_rotated = 0;
   boosted250_decay_id = 0;
   boosted250_decay_isstar = 0;
   boosted250_decay_h_dr = 0;
   boosted250_decay_h_deta = 0;
   boosted250_decay_h_dphi = 0;
   boosted250_decay_h_deta_rotated = 0;
   boosted250_decay_h_dphi_rotated = 0;
   boosted250_lepton_id = 0;
   boosted250_lepton_isstar = 0;
   boosted250_lepton_h_dr = 0;
   boosted250_lepton_h_deta = 0;
   boosted250_lepton_h_dphi = 0;
   boosted250_lepton_h_deta_rotated = 0;
   boosted250_lepton_h_dphi_rotated = 0;
   boosted250_quark_id = 0;
   boosted250_quark_isstar = 0;
   boosted250_quark_h_dr = 0;
   boosted250_quark_h_deta = 0;
   boosted250_quark_h_dphi = 0;
   boosted250_quark_h_deta_rotated = 0;
   boosted250_quark_h_dphi_rotated = 0;
   boosted500_decay_id = 0;
   boosted500_decay_isstar = 0;
   boosted500_decay_h_dr = 0;
   boosted500_decay_h_deta = 0;
   boosted500_decay_h_dphi = 0;
   boosted500_decay_h_deta_rotated = 0;
   boosted500_decay_h_dphi_rotated = 0;
   boosted500_lepton_id = 0;
   boosted500_lepton_isstar = 0;
   boosted500_lepton_h_dr = 0;
   boosted500_lepton_h_deta = 0;
   boosted500_lepton_h_dphi = 0;
   boosted500_lepton_h_deta_rotated = 0;
   boosted500_lepton_h_dphi_rotated = 0;
   boosted500_quark_id = 0;
   boosted500_quark_isstar = 0;
   boosted500_quark_h_dr = 0;
   boosted500_quark_h_deta = 0;
   boosted500_quark_h_dphi = 0;
   boosted500_quark_h_deta_rotated = 0;
   boosted500_quark_h_dphi_rotated = 0;
   boosted1000_decay_id = 0;
   boosted1000_decay_isstar = 0;
   boosted1000_decay_h_dr = 0;
   boosted1000_decay_h_deta = 0;
   boosted1000_decay_h_dphi = 0;
   boosted1000_decay_h_deta_rotated = 0;
   boosted1000_decay_h_dphi_rotated = 0;
   boosted1000_lepton_id = 0;
   boosted1000_lepton_isstar = 0;
   boosted1000_lepton_h_dr = 0;
   boosted1000_lepton_h_deta = 0;
   boosted1000_lepton_h_dphi = 0;
   boosted1000_lepton_h_deta_rotated = 0;
   boosted1000_lepton_h_dphi_rotated = 0;
   boosted1000_quark_id = 0;
   boosted1000_quark_isstar = 0;
   boosted1000_quark_h_dr = 0;
   boosted1000_quark_h_deta = 0;
   boosted1000_quark_h_dphi = 0;
   boosted1000_quark_h_deta_rotated = 0;
   boosted1000_quark_h_dphi_rotated = 0;
   boosted1500_decay_id = 0;
   boosted1500_decay_isstar = 0;
   boosted1500_decay_h_dr = 0;
   boosted1500_decay_h_deta = 0;
   boosted1500_decay_h_dphi = 0;
   boosted1500_decay_h_deta_rotated = 0;
   boosted1500_decay_h_dphi_rotated = 0;
   boosted1500_lepton_id = 0;
   boosted1500_lepton_isstar = 0;
   boosted1500_lepton_h_dr = 0;
   boosted1500_lepton_h_deta = 0;
   boosted1500_lepton_h_dphi = 0;
   boosted1500_lepton_h_deta_rotated = 0;
   boosted1500_lepton_h_dphi_rotated = 0;
   boosted1500_quark_id = 0;
   boosted1500_quark_isstar = 0;
   boosted1500_quark_h_dr = 0;
   boosted1500_quark_h_deta = 0;
   boosted1500_quark_h_dphi = 0;
   boosted1500_quark_h_deta_rotated = 0;
   boosted1500_quark_h_dphi_rotated = 0;
   w_islep = 0;
   w_isstar = 0;
   w_isH = 0;
   l_w_pt = 0;
   l_w_eta = 0;
   l_w_phi = 0;
   l_w_mass = 0;
   l_w_id = 0;
   l_isstar = 0;
   l_isH = 0;
   l_istau = 0;
   q_w_pt = 0;
   q_w_eta = 0;
   q_w_phi = 0;
   q_w_mass = 0;
   q_w_id = 0;
   q_isstar = 0;
   q_isH = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("evt_scale1fb", &evt_scale1fb, &b_evt_scale1fb);
   fChain->SetBranchAddress("xsec_br", &xsec_br, &b_xsec_br);
   fChain->SetBranchAddress("evt_passgoodrunlist", &evt_passgoodrunlist, &b_evt_passgoodrunlist);
   fChain->SetBranchAddress("CMS4path", &CMS4path, &b_CMS4path);
   fChain->SetBranchAddress("CMS4index", &CMS4index, &b_CMS4index);
   fChain->SetBranchAddress("weight_fr_r1_f1", &weight_fr_r1_f1, &b_weight_fr_r1_f1);
   fChain->SetBranchAddress("weight_fr_r1_f2", &weight_fr_r1_f2, &b_weight_fr_r1_f2);
   fChain->SetBranchAddress("weight_fr_r1_f0p5", &weight_fr_r1_f0p5, &b_weight_fr_r1_f0p5);
   fChain->SetBranchAddress("weight_fr_r2_f1", &weight_fr_r2_f1, &b_weight_fr_r2_f1);
   fChain->SetBranchAddress("weight_fr_r2_f2", &weight_fr_r2_f2, &b_weight_fr_r2_f2);
   fChain->SetBranchAddress("weight_fr_r2_f0p5", &weight_fr_r2_f0p5, &b_weight_fr_r2_f0p5);
   fChain->SetBranchAddress("weight_fr_r0p5_f1", &weight_fr_r0p5_f1, &b_weight_fr_r0p5_f1);
   fChain->SetBranchAddress("weight_fr_r0p5_f2", &weight_fr_r0p5_f2, &b_weight_fr_r0p5_f2);
   fChain->SetBranchAddress("weight_fr_r0p5_f0p5", &weight_fr_r0p5_f0p5, &b_weight_fr_r0p5_f0p5);
   fChain->SetBranchAddress("weight_pdf_up", &weight_pdf_up, &b_weight_pdf_up);
   fChain->SetBranchAddress("weight_pdf_down", &weight_pdf_down, &b_weight_pdf_down);
   fChain->SetBranchAddress("weight_alphas_down", &weight_alphas_down, &b_weight_alphas_down);
   fChain->SetBranchAddress("weight_alphas_up", &weight_alphas_up, &b_weight_alphas_up);
   fChain->SetBranchAddress("weight_isr", &weight_isr, &b_weight_isr);
   fChain->SetBranchAddress("weight_isr_up", &weight_isr_up, &b_weight_isr_up);
   fChain->SetBranchAddress("weight_isr_down", &weight_isr_down, &b_weight_isr_down);
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_DoubleEl_DZ", &HLT_DoubleEl_DZ, &b_HLT_DoubleEl_DZ);
   fChain->SetBranchAddress("HLT_DoubleEl_DZ_2", &HLT_DoubleEl_DZ_2, &b_HLT_DoubleEl_DZ_2);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("HLT_SingleEl8", &HLT_SingleEl8, &b_HLT_SingleEl8);
   fChain->SetBranchAddress("HLT_SingleEl17", &HLT_SingleEl17, &b_HLT_SingleEl17);
   fChain->SetBranchAddress("HLT_SingleIsoEl8", &HLT_SingleIsoEl8, &b_HLT_SingleIsoEl8);
   fChain->SetBranchAddress("HLT_SingleIsoEl17", &HLT_SingleIsoEl17, &b_HLT_SingleIsoEl17);
   fChain->SetBranchAddress("HLT_SingleIsoEl23", &HLT_SingleIsoEl23, &b_HLT_SingleIsoEl23);
   fChain->SetBranchAddress("HLT_SingleIsoMu8", &HLT_SingleIsoMu8, &b_HLT_SingleIsoMu8);
   fChain->SetBranchAddress("HLT_SingleIsoMu17", &HLT_SingleIsoMu17, &b_HLT_SingleIsoMu17);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("mc_HLT_DoubleMu", &mc_HLT_DoubleMu, &b_mc_HLT_DoubleMu);
   fChain->SetBranchAddress("mc_HLT_DoubleEl", &mc_HLT_DoubleEl, &b_mc_HLT_DoubleEl);
   fChain->SetBranchAddress("mc_HLT_DoubleEl_DZ", &mc_HLT_DoubleEl_DZ, &b_mc_HLT_DoubleEl_DZ);
   fChain->SetBranchAddress("mc_HLT_DoubleEl_DZ_2", &mc_HLT_DoubleEl_DZ_2, &b_mc_HLT_DoubleEl_DZ_2);
   fChain->SetBranchAddress("mc_HLT_MuEG", &mc_HLT_MuEG, &b_mc_HLT_MuEG);
   fChain->SetBranchAddress("mc_HLT_SingleEl8", &mc_HLT_SingleEl8, &b_mc_HLT_SingleEl8);
   fChain->SetBranchAddress("mc_HLT_SingleEl17", &mc_HLT_SingleEl17, &b_mc_HLT_SingleEl17);
   fChain->SetBranchAddress("mc_HLT_SingleIsoEl8", &mc_HLT_SingleIsoEl8, &b_mc_HLT_SingleIsoEl8);
   fChain->SetBranchAddress("mc_HLT_SingleIsoEl17", &mc_HLT_SingleIsoEl17, &b_mc_HLT_SingleIsoEl17);
   fChain->SetBranchAddress("mc_HLT_SingleIsoEl23", &mc_HLT_SingleIsoEl23, &b_mc_HLT_SingleIsoEl23);
   fChain->SetBranchAddress("mc_HLT_SingleIsoMu8", &mc_HLT_SingleIsoMu8, &b_mc_HLT_SingleIsoMu8);
   fChain->SetBranchAddress("mc_HLT_SingleIsoMu17", &mc_HLT_SingleIsoMu17, &b_mc_HLT_SingleIsoMu17);
   fChain->SetBranchAddress("mc_HLT_PFMET140_PFMHT140_IDTight", &mc_HLT_PFMET140_PFMHT140_IDTight, &b_mc_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("pass_duplicate_ee_em_mm", &pass_duplicate_ee_em_mm, &b_pass_duplicate_ee_em_mm);
   fChain->SetBranchAddress("pass_duplicate_mm_em_ee", &pass_duplicate_mm_em_ee, &b_pass_duplicate_mm_em_ee);
   fChain->SetBranchAddress("is2016", &is2016, &b_is2016);
   fChain->SetBranchAddress("is2017", &is2017, &b_is2017);
   fChain->SetBranchAddress("HLT_MuEG_2016", &HLT_MuEG_2016, &b_HLT_MuEG_2016);
   fChain->SetBranchAddress("mc_HLT_MuEG_2016", &mc_HLT_MuEG_2016, &b_mc_HLT_MuEG_2016);
   fChain->SetBranchAddress("pass_duplicate_ee_em2016_mm", &pass_duplicate_ee_em2016_mm, &b_pass_duplicate_ee_em2016_mm);
   fChain->SetBranchAddress("pass_duplicate_mm_em2016_ee", &pass_duplicate_mm_em2016_ee, &b_pass_duplicate_mm_em2016_ee);
   fChain->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger);
   fChain->SetBranchAddress("lep_p4", &lep_p4_, &b_lep_p4_);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fX", lep_p4_fCoordinates_fX, &b_lep_p4_fCoordinates_fX);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fY", lep_p4_fCoordinates_fY, &b_lep_p4_fCoordinates_fY);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fZ", lep_p4_fCoordinates_fZ, &b_lep_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fT", lep_p4_fCoordinates_fT, &b_lep_p4_fCoordinates_fT);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_coneCorrPt", &lep_coneCorrPt, &b_lep_coneCorrPt);
   fChain->SetBranchAddress("lep_ip3d", &lep_ip3d, &b_lep_ip3d);
   fChain->SetBranchAddress("lep_ip3derr", &lep_ip3derr, &b_lep_ip3derr);
   fChain->SetBranchAddress("lep_isTriggerSafe_v1", &lep_isTriggerSafe_v1, &b_lep_isTriggerSafe_v1);
   fChain->SetBranchAddress("lep_lostHits", &lep_lostHits, &b_lep_lostHits);
   fChain->SetBranchAddress("lep_convVeto", &lep_convVeto, &b_lep_convVeto);
   fChain->SetBranchAddress("lep_motherIdSS", &lep_motherIdSS, &b_lep_motherIdSS);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_3l_fo", &lep_pass_VVV_cutbased_3l_fo, &b_lep_pass_VVV_cutbased_3l_fo);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_3l_tight", &lep_pass_VVV_cutbased_3l_tight, &b_lep_pass_VVV_cutbased_3l_tight);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_fo", &lep_pass_VVV_cutbased_fo, &b_lep_pass_VVV_cutbased_fo);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_tight", &lep_pass_VVV_cutbased_tight, &b_lep_pass_VVV_cutbased_tight);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_veto", &lep_pass_VVV_cutbased_veto, &b_lep_pass_VVV_cutbased_veto);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_fo_noiso", &lep_pass_VVV_cutbased_fo_noiso, &b_lep_pass_VVV_cutbased_fo_noiso);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_tight_noiso", &lep_pass_VVV_cutbased_tight_noiso, &b_lep_pass_VVV_cutbased_tight_noiso);
   fChain->SetBranchAddress("lep_pass_VVV_cutbased_veto_noiso", &lep_pass_VVV_cutbased_veto_noiso, &b_lep_pass_VVV_cutbased_veto_noiso);
   fChain->SetBranchAddress("lep_pass_POG_veto", &lep_pass_POG_veto, &b_lep_pass_POG_veto);
   fChain->SetBranchAddress("lep_pass_POG_loose", &lep_pass_POG_loose, &b_lep_pass_POG_loose);
   fChain->SetBranchAddress("lep_pass_POG_medium", &lep_pass_POG_medium, &b_lep_pass_POG_medium);
   fChain->SetBranchAddress("lep_pass_POG_tight", &lep_pass_POG_tight, &b_lep_pass_POG_tight);
   fChain->SetBranchAddress("lep_pdgId", &lep_pdgId, &b_lep_pdgId);
   fChain->SetBranchAddress("lep_dxy", &lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", &lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_pterr", &lep_pterr, &b_lep_pterr);
   fChain->SetBranchAddress("lep_relIso04DB", &lep_relIso04DB, &b_lep_relIso04DB);
   fChain->SetBranchAddress("lep_relIso03EA", &lep_relIso03EA, &b_lep_relIso03EA);
   fChain->SetBranchAddress("lep_relIso03EALep", &lep_relIso03EALep, &b_lep_relIso03EALep);
   fChain->SetBranchAddress("lep_relIso03EAv2", &lep_relIso03EAv2, &b_lep_relIso03EAv2);
   fChain->SetBranchAddress("lep_relIso04EAv2", &lep_relIso04EAv2, &b_lep_relIso04EAv2);
   fChain->SetBranchAddress("lep_relIso03EAv2Lep", &lep_relIso03EAv2Lep, &b_lep_relIso03EAv2Lep);
   fChain->SetBranchAddress("lep_tightCharge", &lep_tightCharge, &b_lep_tightCharge);
   fChain->SetBranchAddress("lep_trk_pt", &lep_trk_pt, &b_lep_trk_pt);
   fChain->SetBranchAddress("lep_charge", &lep_charge, &b_lep_charge);
   fChain->SetBranchAddress("lep_etaSC", &lep_etaSC, &b_lep_etaSC);
   fChain->SetBranchAddress("lep_MVA", &lep_MVA, &b_lep_MVA);
   fChain->SetBranchAddress("lep_isMediumPOG", &lep_isMediumPOG, &b_lep_isMediumPOG);
   fChain->SetBranchAddress("lep_isTightPOG", &lep_isTightPOG, &b_lep_isTightPOG);
   fChain->SetBranchAddress("lep_isFromW", &lep_isFromW, &b_lep_isFromW);
   fChain->SetBranchAddress("lep_isFromZ", &lep_isFromZ, &b_lep_isFromZ);
   fChain->SetBranchAddress("lep_isFromB", &lep_isFromB, &b_lep_isFromB);
   fChain->SetBranchAddress("lep_isFromC", &lep_isFromC, &b_lep_isFromC);
   fChain->SetBranchAddress("lep_isFromL", &lep_isFromL, &b_lep_isFromL);
   fChain->SetBranchAddress("lep_isFromLF", &lep_isFromLF, &b_lep_isFromLF);
   fChain->SetBranchAddress("lep_genPart_index", &lep_genPart_index, &b_lep_genPart_index);
   fChain->SetBranchAddress("lep_r9", &lep_r9, &b_lep_r9);
   fChain->SetBranchAddress("lep_nlayers", &lep_nlayers, &b_lep_nlayers);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_relIso03EA", &el_relIso03EA, &b_el_relIso03EA);
   fChain->SetBranchAddress("el_relIso03EALep", &el_relIso03EALep, &b_el_relIso03EALep);
   fChain->SetBranchAddress("el_ip3d", &el_ip3d, &b_el_ip3d);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_relIso04DB", &mu_relIso04DB, &b_mu_relIso04DB);
   fChain->SetBranchAddress("mu_relIso03EA", &mu_relIso03EA, &b_mu_relIso03EA);
   fChain->SetBranchAddress("mu_relIso03EALep", &mu_relIso03EALep, &b_mu_relIso03EALep);
   fChain->SetBranchAddress("mu_ip3d", &mu_ip3d, &b_mu_ip3d);
   fChain->SetBranchAddress("lbnt_pt", &lbnt_pt, &b_lbnt_pt);
   fChain->SetBranchAddress("lbnt_coneCorrPt", &lbnt_coneCorrPt, &b_lbnt_coneCorrPt);
   fChain->SetBranchAddress("lbnt_abseta", &lbnt_abseta, &b_lbnt_abseta);
   fChain->SetBranchAddress("lbnt_pdgId", &lbnt_pdgId, &b_lbnt_pdgId);
   fChain->SetBranchAddress("lbnt_el_pt", &lbnt_el_pt, &b_lbnt_el_pt);
   fChain->SetBranchAddress("lbnt_el_coneCorrPt", &lbnt_el_coneCorrPt, &b_lbnt_el_coneCorrPt);
   fChain->SetBranchAddress("lbnt_el_abseta", &lbnt_el_abseta, &b_lbnt_el_abseta);
   fChain->SetBranchAddress("lbnt_mu_pt", &lbnt_mu_pt, &b_lbnt_mu_pt);
   fChain->SetBranchAddress("lbnt_mu_coneCorrPt", &lbnt_mu_coneCorrPt, &b_lbnt_mu_coneCorrPt);
   fChain->SetBranchAddress("lbnt_mu_abseta", &lbnt_mu_abseta, &b_lbnt_mu_abseta);
   fChain->SetBranchAddress("jets_p4", &jets_p4_, &b_jets_p4_);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fX", jets_p4_fCoordinates_fX, &b_jets_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fY", jets_p4_fCoordinates_fY, &b_jets_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fZ", jets_p4_fCoordinates_fZ, &b_jets_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fT", jets_p4_fCoordinates_fT, &b_jets_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_up_p4", &jets_up_p4_, &b_jets_up_p4_);
   fChain->SetBranchAddress("jets_up_p4.fCoordinates.fX", jets_up_p4_fCoordinates_fX, &b_jets_up_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_up_p4.fCoordinates.fY", jets_up_p4_fCoordinates_fY, &b_jets_up_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_up_p4.fCoordinates.fZ", jets_up_p4_fCoordinates_fZ, &b_jets_up_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_up_p4.fCoordinates.fT", jets_up_p4_fCoordinates_fT, &b_jets_up_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_dn_p4", &jets_dn_p4_, &b_jets_dn_p4_);
   fChain->SetBranchAddress("jets_dn_p4.fCoordinates.fX", jets_dn_p4_fCoordinates_fX, &b_jets_dn_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_dn_p4.fCoordinates.fY", jets_dn_p4_fCoordinates_fY, &b_jets_dn_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_dn_p4.fCoordinates.fZ", jets_dn_p4_fCoordinates_fZ, &b_jets_dn_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_dn_p4.fCoordinates.fT", jets_dn_p4_fCoordinates_fT, &b_jets_dn_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_csv", &jets_csv, &b_jets_csv);
   fChain->SetBranchAddress("jets_up_csv", &jets_up_csv, &b_jets_up_csv);
   fChain->SetBranchAddress("jets_dn_csv", &jets_dn_csv, &b_jets_dn_csv);
   fChain->SetBranchAddress("jets_jer_csv", &jets_jer_csv, &b_jets_jer_csv);
   fChain->SetBranchAddress("jets_jerup_csv", &jets_jerup_csv, &b_jets_jerup_csv);
   fChain->SetBranchAddress("jets_jerdn_csv", &jets_jerdn_csv, &b_jets_jerdn_csv);
   fChain->SetBranchAddress("jets_jer_p4", &jets_jer_p4_, &b_jets_jer_p4_);
   fChain->SetBranchAddress("jets_jer_p4.fCoordinates.fX", jets_jer_p4_fCoordinates_fX, &b_jets_jer_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_jer_p4.fCoordinates.fY", jets_jer_p4_fCoordinates_fY, &b_jets_jer_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_jer_p4.fCoordinates.fZ", jets_jer_p4_fCoordinates_fZ, &b_jets_jer_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_jer_p4.fCoordinates.fT", jets_jer_p4_fCoordinates_fT, &b_jets_jer_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_jerup_p4", &jets_jerup_p4_, &b_jets_jerup_p4_);
   fChain->SetBranchAddress("jets_jerup_p4.fCoordinates.fX", jets_jerup_p4_fCoordinates_fX, &b_jets_jerup_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_jerup_p4.fCoordinates.fY", jets_jerup_p4_fCoordinates_fY, &b_jets_jerup_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_jerup_p4.fCoordinates.fZ", jets_jerup_p4_fCoordinates_fZ, &b_jets_jerup_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_jerup_p4.fCoordinates.fT", jets_jerup_p4_fCoordinates_fT, &b_jets_jerup_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_jerdn_p4", &jets_jerdn_p4_, &b_jets_jerdn_p4_);
   fChain->SetBranchAddress("jets_jerdn_p4.fCoordinates.fX", jets_jerdn_p4_fCoordinates_fX, &b_jets_jerdn_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_jerdn_p4.fCoordinates.fY", jets_jerdn_p4_fCoordinates_fY, &b_jets_jerdn_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_jerdn_p4.fCoordinates.fZ", jets_jerdn_p4_fCoordinates_fZ, &b_jets_jerdn_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_jerdn_p4.fCoordinates.fT", jets_jerdn_p4_fCoordinates_fT, &b_jets_jerdn_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_p4", &jets30_p4_, &b_jets30_p4_);
   fChain->SetBranchAddress("jets30_p4.fCoordinates.fX", jets30_p4_fCoordinates_fX, &b_jets30_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_p4.fCoordinates.fY", jets30_p4_fCoordinates_fY, &b_jets30_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_p4.fCoordinates.fZ", jets30_p4_fCoordinates_fZ, &b_jets30_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_p4.fCoordinates.fT", jets30_p4_fCoordinates_fT, &b_jets30_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_up_p4", &jets30_up_p4_, &b_jets30_up_p4_);
   fChain->SetBranchAddress("jets30_up_p4.fCoordinates.fX", jets30_up_p4_fCoordinates_fX, &b_jets30_up_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_up_p4.fCoordinates.fY", jets30_up_p4_fCoordinates_fY, &b_jets30_up_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_up_p4.fCoordinates.fZ", jets30_up_p4_fCoordinates_fZ, &b_jets30_up_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_up_p4.fCoordinates.fT", jets30_up_p4_fCoordinates_fT, &b_jets30_up_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_dn_p4", &jets30_dn_p4_, &b_jets30_dn_p4_);
   fChain->SetBranchAddress("jets30_dn_p4.fCoordinates.fX", jets30_dn_p4_fCoordinates_fX, &b_jets30_dn_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_dn_p4.fCoordinates.fY", jets30_dn_p4_fCoordinates_fY, &b_jets30_dn_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_dn_p4.fCoordinates.fZ", jets30_dn_p4_fCoordinates_fZ, &b_jets30_dn_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_dn_p4.fCoordinates.fT", jets30_dn_p4_fCoordinates_fT, &b_jets30_dn_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_jer_p4", &jets30_jer_p4_, &b_jets30_jer_p4_);
   fChain->SetBranchAddress("jets30_jer_p4.fCoordinates.fX", jets30_jer_p4_fCoordinates_fX, &b_jets30_jer_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_jer_p4.fCoordinates.fY", jets30_jer_p4_fCoordinates_fY, &b_jets30_jer_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_jer_p4.fCoordinates.fZ", jets30_jer_p4_fCoordinates_fZ, &b_jets30_jer_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_jer_p4.fCoordinates.fT", jets30_jer_p4_fCoordinates_fT, &b_jets30_jer_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_jerup_p4", &jets30_jerup_p4_, &b_jets30_jerup_p4_);
   fChain->SetBranchAddress("jets30_jerup_p4.fCoordinates.fX", jets30_jerup_p4_fCoordinates_fX, &b_jets30_jerup_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_jerup_p4.fCoordinates.fY", jets30_jerup_p4_fCoordinates_fY, &b_jets30_jerup_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_jerup_p4.fCoordinates.fZ", jets30_jerup_p4_fCoordinates_fZ, &b_jets30_jerup_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_jerup_p4.fCoordinates.fT", jets30_jerup_p4_fCoordinates_fT, &b_jets30_jerup_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets30_jerdn_p4", &jets30_jerdn_p4_, &b_jets30_jerdn_p4_);
   fChain->SetBranchAddress("jets30_jerdn_p4.fCoordinates.fX", jets30_jerdn_p4_fCoordinates_fX, &b_jets30_jerdn_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets30_jerdn_p4.fCoordinates.fY", jets30_jerdn_p4_fCoordinates_fY, &b_jets30_jerdn_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets30_jerdn_p4.fCoordinates.fZ", jets30_jerdn_p4_fCoordinates_fZ, &b_jets30_jerdn_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets30_jerdn_p4.fCoordinates.fT", jets30_jerdn_p4_fCoordinates_fT, &b_jets30_jerdn_p4_fCoordinates_fT);
   fChain->SetBranchAddress("ak8jets_p4", &ak8jets_p4_, &b_ak8jets_p4_);
   fChain->SetBranchAddress("ak8jets_p4.fCoordinates.fX", ak8jets_p4_fCoordinates_fX, &b_ak8jets_p4_fCoordinates_fX);
   fChain->SetBranchAddress("ak8jets_p4.fCoordinates.fY", ak8jets_p4_fCoordinates_fY, &b_ak8jets_p4_fCoordinates_fY);
   fChain->SetBranchAddress("ak8jets_p4.fCoordinates.fZ", ak8jets_p4_fCoordinates_fZ, &b_ak8jets_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("ak8jets_p4.fCoordinates.fT", ak8jets_p4_fCoordinates_fT, &b_ak8jets_p4_fCoordinates_fT);
   fChain->SetBranchAddress("ak8jets_softdropMass", &ak8jets_softdropMass, &b_ak8jets_softdropMass);
   fChain->SetBranchAddress("ak8jets_prunedMass", &ak8jets_prunedMass, &b_ak8jets_prunedMass);
   fChain->SetBranchAddress("ak8jets_trimmedMass", &ak8jets_trimmedMass, &b_ak8jets_trimmedMass);
   fChain->SetBranchAddress("ak8jets_mass", &ak8jets_mass, &b_ak8jets_mass);
   fChain->SetBranchAddress("ak8jets_nJettinessTau1", &ak8jets_nJettinessTau1, &b_ak8jets_nJettinessTau1);
   fChain->SetBranchAddress("ak8jets_nJettinessTau2", &ak8jets_nJettinessTau2, &b_ak8jets_nJettinessTau2);
   fChain->SetBranchAddress("ak8jets_softdropPuppiSubjet1", &ak8jets_softdropPuppiSubjet1, &b_ak8jets_softdropPuppiSubjet1);
   fChain->SetBranchAddress("ak8jets_softdropPuppiSubjet2", &ak8jets_softdropPuppiSubjet2, &b_ak8jets_softdropPuppiSubjet2);
   fChain->SetBranchAddress("ak8jets_puppi_softdropMass", &ak8jets_puppi_softdropMass, &b_ak8jets_puppi_softdropMass);
   fChain->SetBranchAddress("ak8jets_puppi_nJettinessTau1", &ak8jets_puppi_nJettinessTau1, &b_ak8jets_puppi_nJettinessTau1);
   fChain->SetBranchAddress("ak8jets_puppi_nJettinessTau2", &ak8jets_puppi_nJettinessTau2, &b_ak8jets_puppi_nJettinessTau2);
   fChain->SetBranchAddress("ak8jets_puppi_eta", &ak8jets_puppi_eta, &b_ak8jets_puppi_eta);
   fChain->SetBranchAddress("ak8jets_puppi_phi", &ak8jets_puppi_phi, &b_ak8jets_puppi_phi);
   fChain->SetBranchAddress("ak8jets_puppi_pt", &ak8jets_puppi_pt, &b_ak8jets_puppi_pt);
   fChain->SetBranchAddress("ak8jets_puppi_mass", &ak8jets_puppi_mass, &b_ak8jets_puppi_mass);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_up_pt", &met_up_pt, &b_met_up_pt);
   fChain->SetBranchAddress("met_up_phi", &met_up_phi, &b_met_up_phi);
   fChain->SetBranchAddress("met_dn_pt", &met_dn_pt, &b_met_dn_pt);
   fChain->SetBranchAddress("met_dn_phi", &met_dn_phi, &b_met_dn_phi);
   fChain->SetBranchAddress("met_gen_pt", &met_gen_pt, &b_met_gen_pt);
   fChain->SetBranchAddress("met_gen_phi", &met_gen_phi, &b_met_gen_phi);
   fChain->SetBranchAddress("met_jer_pt", &met_jer_pt, &b_met_jer_pt);
   fChain->SetBranchAddress("met_jerup_pt", &met_jerup_pt, &b_met_jerup_pt);
   fChain->SetBranchAddress("met_jerdn_pt", &met_jerdn_pt, &b_met_jerdn_pt);
   fChain->SetBranchAddress("met_jer_phi", &met_jer_phi, &b_met_jer_phi);
   fChain->SetBranchAddress("met_jerup_phi", &met_jerup_phi, &b_met_jerup_phi);
   fChain->SetBranchAddress("met_jerdn_phi", &met_jerdn_phi, &b_met_jerdn_phi);
   fChain->SetBranchAddress("firstgoodvertex", &firstgoodvertex, &b_firstgoodvertex);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("nVert", &nVert, &b_nVert);
   fChain->SetBranchAddress("nisoTrack_mt2_cleaned_VVV_cutbased_veto", &nisoTrack_mt2_cleaned_VVV_cutbased_veto, &b_nisoTrack_mt2_cleaned_VVV_cutbased_veto);
   fChain->SetBranchAddress("weight_btagsf", &weight_btagsf, &b_weight_btagsf);
   fChain->SetBranchAddress("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, &b_weight_btagsf_heavy_DN);
   fChain->SetBranchAddress("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, &b_weight_btagsf_heavy_UP);
   fChain->SetBranchAddress("weight_btagsf_light_DN", &weight_btagsf_light_DN, &b_weight_btagsf_light_DN);
   fChain->SetBranchAddress("weight_btagsf_light_UP", &weight_btagsf_light_UP, &b_weight_btagsf_light_UP);
   fChain->SetBranchAddress("gen_ht", &gen_ht, &b_gen_ht);
   fChain->SetBranchAddress("genPart_p4", &genPart_p4_, &b_genPart_p4_);
   fChain->SetBranchAddress("genPart_p4.fCoordinates.fX", genPart_p4_fCoordinates_fX, &b_genPart_p4_fCoordinates_fX);
   fChain->SetBranchAddress("genPart_p4.fCoordinates.fY", genPart_p4_fCoordinates_fY, &b_genPart_p4_fCoordinates_fY);
   fChain->SetBranchAddress("genPart_p4.fCoordinates.fZ", genPart_p4_fCoordinates_fZ, &b_genPart_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("genPart_p4.fCoordinates.fT", genPart_p4_fCoordinates_fT, &b_genPart_p4_fCoordinates_fT);
   fChain->SetBranchAddress("genPart_motherId", &genPart_motherId, &b_genPart_motherId);
   fChain->SetBranchAddress("genPart_pdgId", &genPart_pdgId, &b_genPart_pdgId);
   fChain->SetBranchAddress("genPart_charge", &genPart_charge, &b_genPart_charge);
   fChain->SetBranchAddress("genPart_status", &genPart_status, &b_genPart_status);
   fChain->SetBranchAddress("ngenLep", &ngenLep, &b_ngenLep);
   fChain->SetBranchAddress("ngenLepFromTau", &ngenLepFromTau, &b_ngenLepFromTau);
   fChain->SetBranchAddress("ngenLepFromBoson", &ngenLepFromBoson, &b_ngenLepFromBoson);
   fChain->SetBranchAddress("Flag_AllEventFilters", &Flag_AllEventFilters, &b_Flag_AllEventFilters);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, &b_Flag_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, &b_Flag_badChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_badMuonFilter", &Flag_badMuonFilter, &b_Flag_badMuonFilter);
   fChain->SetBranchAddress("Flag_badMuonFilterv2", &Flag_badMuonFilterv2, &b_Flag_badMuonFilterv2);
   fChain->SetBranchAddress("Flag_badChargedCandidateFilterv2", &Flag_badChargedCandidateFilterv2, &b_Flag_badChargedCandidateFilterv2);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016", &Flag_globalTightHalo2016, &b_Flag_globalTightHalo2016);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, &b_Flag_trackingFailureFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_badMuons", &Flag_badMuons, &b_Flag_badMuons);
   fChain->SetBranchAddress("Flag_duplicateMuons", &Flag_duplicateMuons, &b_Flag_duplicateMuons);
   fChain->SetBranchAddress("Flag_noBadMuons", &Flag_noBadMuons, &b_Flag_noBadMuons);
   fChain->SetBranchAddress("fastsimfilt", &fastsimfilt, &b_fastsimfilt);
   fChain->SetBranchAddress("nVlep", &nVlep, &b_nVlep);
   fChain->SetBranchAddress("nTlep", &nTlep, &b_nTlep);
   fChain->SetBranchAddress("nTlepSS", &nTlepSS, &b_nTlepSS);
   fChain->SetBranchAddress("nLlep", &nLlep, &b_nLlep);
   fChain->SetBranchAddress("nLlep3L", &nLlep3L, &b_nLlep3L);
   fChain->SetBranchAddress("nTlep3L", &nTlep3L, &b_nTlep3L);
   fChain->SetBranchAddress("nSFOS", &nSFOS, &b_nSFOS);
   fChain->SetBranchAddress("nSFOSinZ", &nSFOSinZ, &b_nSFOSinZ);
   fChain->SetBranchAddress("nj", &nj, &b_nj);
   fChain->SetBranchAddress("nj_up", &nj_up, &b_nj_up);
   fChain->SetBranchAddress("nj_dn", &nj_dn, &b_nj_dn);
   fChain->SetBranchAddress("nj_jer", &nj_jer, &b_nj_jer);
   fChain->SetBranchAddress("nj_jerup", &nj_jerup, &b_nj_jerup);
   fChain->SetBranchAddress("nj_jerdn", &nj_jerdn, &b_nj_jerdn);
   fChain->SetBranchAddress("nj30", &nj30, &b_nj30);
   fChain->SetBranchAddress("nj30_up", &nj30_up, &b_nj30_up);
   fChain->SetBranchAddress("nj30_dn", &nj30_dn, &b_nj30_dn);
   fChain->SetBranchAddress("nj30_jer", &nj30_jer, &b_nj30_jer);
   fChain->SetBranchAddress("nj30_jerup", &nj30_jerup, &b_nj30_jerup);
   fChain->SetBranchAddress("nj30_jerdn", &nj30_jerdn, &b_nj30_jerdn);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("nb_up", &nb_up, &b_nb_up);
   fChain->SetBranchAddress("nb_dn", &nb_dn, &b_nb_dn);
   fChain->SetBranchAddress("nb_jer", &nb_jer, &b_nb_jer);
   fChain->SetBranchAddress("nb_jerup", &nb_jerup, &b_nb_jerup);
   fChain->SetBranchAddress("nb_jerdn", &nb_jerdn, &b_nb_jerdn);
   fChain->SetBranchAddress("Ml0j0", &Ml0j0, &b_Ml0j0);
   fChain->SetBranchAddress("Ml0j0_up", &Ml0j0_up, &b_Ml0j0_up);
   fChain->SetBranchAddress("Ml0j0_dn", &Ml0j0_dn, &b_Ml0j0_dn);
   fChain->SetBranchAddress("Ml0j0_jer", &Ml0j0_jer, &b_Ml0j0_jer);
   fChain->SetBranchAddress("Ml0j0_jerup", &Ml0j0_jerup, &b_Ml0j0_jerup);
   fChain->SetBranchAddress("Ml0j0_jerdn", &Ml0j0_jerdn, &b_Ml0j0_jerdn);
   fChain->SetBranchAddress("Ml0j1", &Ml0j1, &b_Ml0j1);
   fChain->SetBranchAddress("Ml0j1_up", &Ml0j1_up, &b_Ml0j1_up);
   fChain->SetBranchAddress("Ml0j1_dn", &Ml0j1_dn, &b_Ml0j1_dn);
   fChain->SetBranchAddress("Ml0j1_jer", &Ml0j1_jer, &b_Ml0j1_jer);
   fChain->SetBranchAddress("Ml0j1_jerup", &Ml0j1_jerup, &b_Ml0j1_jerup);
   fChain->SetBranchAddress("Ml0j1_jerdn", &Ml0j1_jerdn, &b_Ml0j1_jerdn);
   fChain->SetBranchAddress("Ml1j0", &Ml1j0, &b_Ml1j0);
   fChain->SetBranchAddress("Ml1j0_up", &Ml1j0_up, &b_Ml1j0_up);
   fChain->SetBranchAddress("Ml1j0_dn", &Ml1j0_dn, &b_Ml1j0_dn);
   fChain->SetBranchAddress("Ml1j0_jer", &Ml1j0_jer, &b_Ml1j0_jer);
   fChain->SetBranchAddress("Ml1j0_jerup", &Ml1j0_jerup, &b_Ml1j0_jerup);
   fChain->SetBranchAddress("Ml1j0_jerdn", &Ml1j0_jerdn, &b_Ml1j0_jerdn);
   fChain->SetBranchAddress("Ml1j1", &Ml1j1, &b_Ml1j1);
   fChain->SetBranchAddress("Ml1j1_up", &Ml1j1_up, &b_Ml1j1_up);
   fChain->SetBranchAddress("Ml1j1_dn", &Ml1j1_dn, &b_Ml1j1_dn);
   fChain->SetBranchAddress("Ml1j1_jer", &Ml1j1_jer, &b_Ml1j1_jer);
   fChain->SetBranchAddress("Ml1j1_jerup", &Ml1j1_jerup, &b_Ml1j1_jerup);
   fChain->SetBranchAddress("Ml1j1_jerdn", &Ml1j1_jerdn, &b_Ml1j1_jerdn);
   fChain->SetBranchAddress("MinMlj", &MinMlj, &b_MinMlj);
   fChain->SetBranchAddress("MinMlj_up", &MinMlj_up, &b_MinMlj_up);
   fChain->SetBranchAddress("MinMlj_dn", &MinMlj_dn, &b_MinMlj_dn);
   fChain->SetBranchAddress("MinMlj_jer", &MinMlj_jer, &b_MinMlj_jer);
   fChain->SetBranchAddress("MinMlj_jerup", &MinMlj_jerup, &b_MinMlj_jerup);
   fChain->SetBranchAddress("MinMlj_jerdn", &MinMlj_jerdn, &b_MinMlj_jerdn);
   fChain->SetBranchAddress("SumMinMlj01", &SumMinMlj01, &b_SumMinMlj01);
   fChain->SetBranchAddress("SumMinMlj01_up", &SumMinMlj01_up, &b_SumMinMlj01_up);
   fChain->SetBranchAddress("SumMinMlj01_dn", &SumMinMlj01_dn, &b_SumMinMlj01_dn);
   fChain->SetBranchAddress("SumMinMlj01_jer", &SumMinMlj01_jer, &b_SumMinMlj01_jer);
   fChain->SetBranchAddress("SumMinMlj01_jerup", &SumMinMlj01_jerup, &b_SumMinMlj01_jerup);
   fChain->SetBranchAddress("SumMinMlj01_jerdn", &SumMinMlj01_jerdn, &b_SumMinMlj01_jerdn);
   fChain->SetBranchAddress("MaxMlj", &MaxMlj, &b_MaxMlj);
   fChain->SetBranchAddress("MaxMlj_up", &MaxMlj_up, &b_MaxMlj_up);
   fChain->SetBranchAddress("MaxMlj_dn", &MaxMlj_dn, &b_MaxMlj_dn);
   fChain->SetBranchAddress("MaxMlj_jer", &MaxMlj_jer, &b_MaxMlj_jer);
   fChain->SetBranchAddress("MaxMlj_jerup", &MaxMlj_jerup, &b_MaxMlj_jerup);
   fChain->SetBranchAddress("MaxMlj_jerdn", &MaxMlj_jerdn, &b_MaxMlj_jerdn);
   fChain->SetBranchAddress("SumMlj", &SumMlj, &b_SumMlj);
   fChain->SetBranchAddress("SumMlj_up", &SumMlj_up, &b_SumMlj_up);
   fChain->SetBranchAddress("SumMlj_dn", &SumMlj_dn, &b_SumMlj_dn);
   fChain->SetBranchAddress("SumMlj_jer", &SumMlj_jer, &b_SumMlj_jer);
   fChain->SetBranchAddress("SumMlj_jerup", &SumMlj_jerup, &b_SumMlj_jerup);
   fChain->SetBranchAddress("SumMlj_jerdn", &SumMlj_jerdn, &b_SumMlj_jerdn);
   fChain->SetBranchAddress("Ml0jj", &Ml0jj, &b_Ml0jj);
   fChain->SetBranchAddress("Ml0jj_up", &Ml0jj_up, &b_Ml0jj_up);
   fChain->SetBranchAddress("Ml0jj_dn", &Ml0jj_dn, &b_Ml0jj_dn);
   fChain->SetBranchAddress("Ml0jj_jer", &Ml0jj_jer, &b_Ml0jj_jer);
   fChain->SetBranchAddress("Ml0jj_jerup", &Ml0jj_jerup, &b_Ml0jj_jerup);
   fChain->SetBranchAddress("Ml0jj_jerdn", &Ml0jj_jerdn, &b_Ml0jj_jerdn);
   fChain->SetBranchAddress("Ml1jj", &Ml1jj, &b_Ml1jj);
   fChain->SetBranchAddress("Ml1jj_up", &Ml1jj_up, &b_Ml1jj_up);
   fChain->SetBranchAddress("Ml1jj_dn", &Ml1jj_dn, &b_Ml1jj_dn);
   fChain->SetBranchAddress("Ml1jj_jer", &Ml1jj_jer, &b_Ml1jj_jer);
   fChain->SetBranchAddress("Ml1jj_jerup", &Ml1jj_jerup, &b_Ml1jj_jerup);
   fChain->SetBranchAddress("Ml1jj_jerdn", &Ml1jj_jerdn, &b_Ml1jj_jerdn);
   fChain->SetBranchAddress("MinMljj", &MinMljj, &b_MinMljj);
   fChain->SetBranchAddress("MinMljj_up", &MinMljj_up, &b_MinMljj_up);
   fChain->SetBranchAddress("MinMljj_dn", &MinMljj_dn, &b_MinMljj_dn);
   fChain->SetBranchAddress("MinMljj_jer", &MinMljj_jer, &b_MinMljj_jer);
   fChain->SetBranchAddress("MinMljj_jerup", &MinMljj_jerup, &b_MinMljj_jerup);
   fChain->SetBranchAddress("MinMljj_jerdn", &MinMljj_jerdn, &b_MinMljj_jerdn);
   fChain->SetBranchAddress("MaxMljj", &MaxMljj, &b_MaxMljj);
   fChain->SetBranchAddress("MaxMljj_up", &MaxMljj_up, &b_MaxMljj_up);
   fChain->SetBranchAddress("MaxMljj_dn", &MaxMljj_dn, &b_MaxMljj_dn);
   fChain->SetBranchAddress("MaxMljj_jer", &MaxMljj_jer, &b_MaxMljj_jer);
   fChain->SetBranchAddress("MaxMljj_jerup", &MaxMljj_jerup, &b_MaxMljj_jerup);
   fChain->SetBranchAddress("MaxMljj_jerdn", &MaxMljj_jerdn, &b_MaxMljj_jerdn);
   fChain->SetBranchAddress("SumMljj", &SumMljj, &b_SumMljj);
   fChain->SetBranchAddress("SumMljj_up", &SumMljj_up, &b_SumMljj_up);
   fChain->SetBranchAddress("SumMljj_dn", &SumMljj_dn, &b_SumMljj_dn);
   fChain->SetBranchAddress("SumMljj_jer", &SumMljj_jer, &b_SumMljj_jer);
   fChain->SetBranchAddress("SumMljj_jerup", &SumMljj_jerup, &b_SumMljj_jerup);
   fChain->SetBranchAddress("SumMljj_jerdn", &SumMljj_jerdn, &b_SumMljj_jerdn);
   fChain->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   fChain->SetBranchAddress("Mjj_up", &Mjj_up, &b_Mjj_up);
   fChain->SetBranchAddress("Mjj_dn", &Mjj_dn, &b_Mjj_dn);
   fChain->SetBranchAddress("Mjj_jer", &Mjj_jer, &b_Mjj_jer);
   fChain->SetBranchAddress("Mjj_jerup", &Mjj_jerup, &b_Mjj_jerup);
   fChain->SetBranchAddress("Mjj_jerdn", &Mjj_jerdn, &b_Mjj_jerdn);
   fChain->SetBranchAddress("DRjj", &DRjj, &b_DRjj);
   fChain->SetBranchAddress("DRjj_up", &DRjj_up, &b_DRjj_up);
   fChain->SetBranchAddress("DRjj_dn", &DRjj_dn, &b_DRjj_dn);
   fChain->SetBranchAddress("DRjj_jer", &DRjj_jer, &b_DRjj_jer);
   fChain->SetBranchAddress("DRjj_jerup", &DRjj_jerup, &b_DRjj_jerup);
   fChain->SetBranchAddress("DRjj_jerdn", &DRjj_jerdn, &b_DRjj_jerdn);
   fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_fCoordinates_fX);
   fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_fCoordinates_fY);
   fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_up_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_up_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_up_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_up_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_dn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_dn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_dn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_dn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_jer_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_jer_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_jer_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_jer_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_jerup_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_jerup_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_jerup_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_jerup_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_jerdn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_jerdn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_jerdn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_jerdn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_up_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_up_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_up_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_up_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_dn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_dn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_dn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_dn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_jer_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_jer_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_jer_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_jer_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_jerup_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_jerup_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_jerup_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_jerup_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_jerdn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_jerdn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_jerdn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_jerdn_fCoordinates_fT);
   fChain->SetBranchAddress("MjjDR1", &MjjDR1, &b_MjjDR1);
   fChain->SetBranchAddress("MjjDR1_up", &MjjDR1_up, &b_MjjDR1_up);
   fChain->SetBranchAddress("MjjDR1_dn", &MjjDR1_dn, &b_MjjDR1_dn);
   fChain->SetBranchAddress("MjjDR1_jer", &MjjDR1_jer, &b_MjjDR1_jer);
   fChain->SetBranchAddress("MjjDR1_jerup", &MjjDR1_jerup, &b_MjjDR1_jerup);
   fChain->SetBranchAddress("MjjDR1_jerdn", &MjjDR1_jerdn, &b_MjjDR1_jerdn);
   fChain->SetBranchAddress("DRjjDR1", &DRjjDR1, &b_DRjjDR1);
   fChain->SetBranchAddress("DRjjDR1_up", &DRjjDR1_up, &b_DRjjDR1_up);
   fChain->SetBranchAddress("DRjjDR1_dn", &DRjjDR1_dn, &b_DRjjDR1_dn);
   fChain->SetBranchAddress("DRjjDR1_jer", &DRjjDR1_jer, &b_DRjjDR1_jer);
   fChain->SetBranchAddress("DRjjDR1_jerup", &DRjjDR1_jerup, &b_DRjjDR1_jerup);
   fChain->SetBranchAddress("DRjjDR1_jerdn", &DRjjDR1_jerdn, &b_DRjjDR1_jerdn);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_up_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_up_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_up_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_up_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_dn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_dn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_dn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_dn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_jer_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_jer_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_jer_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_jer_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_jerup_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_jerup_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_jerup_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_jerup_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet0_wtag_p4_DR1_jerdn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_up_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_up_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_up_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_up_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_dn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_dn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_dn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_dn_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_jer_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_jer_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_jer_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_jer_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_jerup_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_jerup_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_jerup_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_jerup_fCoordinates_fT);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_jet1_wtag_p4_DR1_jerdn_fCoordinates_fT);
   fChain->SetBranchAddress("MjjVBF", &MjjVBF, &b_MjjVBF);
   fChain->SetBranchAddress("MjjVBF_up", &MjjVBF_up, &b_MjjVBF_up);
   fChain->SetBranchAddress("MjjVBF_dn", &MjjVBF_dn, &b_MjjVBF_dn);
   fChain->SetBranchAddress("MjjVBF_jer", &MjjVBF_jer, &b_MjjVBF_jer);
   fChain->SetBranchAddress("MjjVBF_jerup", &MjjVBF_jerup, &b_MjjVBF_jerup);
   fChain->SetBranchAddress("MjjVBF_jerdn", &MjjVBF_jerdn, &b_MjjVBF_jerdn);
   fChain->SetBranchAddress("DetajjVBF", &DetajjVBF, &b_DetajjVBF);
   fChain->SetBranchAddress("DetajjVBF_up", &DetajjVBF_up, &b_DetajjVBF_up);
   fChain->SetBranchAddress("DetajjVBF_dn", &DetajjVBF_dn, &b_DetajjVBF_dn);
   fChain->SetBranchAddress("DetajjVBF_jer", &DetajjVBF_jer, &b_DetajjVBF_jer);
   fChain->SetBranchAddress("DetajjVBF_jerup", &DetajjVBF_jerup, &b_DetajjVBF_jerup);
   fChain->SetBranchAddress("DetajjVBF_jerdn", &DetajjVBF_jerdn, &b_DetajjVBF_jerdn);
   fChain->SetBranchAddress("MjjL", &MjjL, &b_MjjL);
   fChain->SetBranchAddress("MjjL_up", &MjjL_up, &b_MjjL_up);
   fChain->SetBranchAddress("MjjL_dn", &MjjL_dn, &b_MjjL_dn);
   fChain->SetBranchAddress("MjjL_jer", &MjjL_jer, &b_MjjL_jer);
   fChain->SetBranchAddress("MjjL_jerup", &MjjL_jerup, &b_MjjL_jerup);
   fChain->SetBranchAddress("MjjL_jerdn", &MjjL_jerdn, &b_MjjL_jerdn);
   fChain->SetBranchAddress("DetajjL", &DetajjL, &b_DetajjL);
   fChain->SetBranchAddress("DetajjL_up", &DetajjL_up, &b_DetajjL_up);
   fChain->SetBranchAddress("DetajjL_dn", &DetajjL_dn, &b_DetajjL_dn);
   fChain->SetBranchAddress("DetajjL_jer", &DetajjL_jer, &b_DetajjL_jer);
   fChain->SetBranchAddress("DetajjL_jerup", &DetajjL_jerup, &b_DetajjL_jerup);
   fChain->SetBranchAddress("DetajjL_jerdn", &DetajjL_jerdn, &b_DetajjL_jerdn);
   fChain->SetBranchAddress("MllSS", &MllSS, &b_MllSS);
   fChain->SetBranchAddress("MeeSS", &MeeSS, &b_MeeSS);
   fChain->SetBranchAddress("Mll3L", &Mll3L, &b_Mll3L);
   fChain->SetBranchAddress("Mee3L", &Mee3L, &b_Mee3L);
   fChain->SetBranchAddress("Mll3L1", &Mll3L1, &b_Mll3L1);
   fChain->SetBranchAddress("M3l", &M3l, &b_M3l);
   fChain->SetBranchAddress("Pt3l", &Pt3l, &b_Pt3l);
   fChain->SetBranchAddress("M01", &M01, &b_M01);
   fChain->SetBranchAddress("M02", &M02, &b_M02);
   fChain->SetBranchAddress("M12", &M12, &b_M12);
   fChain->SetBranchAddress("isSFOS01", &isSFOS01, &b_isSFOS01);
   fChain->SetBranchAddress("isSFOS02", &isSFOS02, &b_isSFOS02);
   fChain->SetBranchAddress("isSFOS12", &isSFOS12, &b_isSFOS12);
   fChain->SetBranchAddress("DPhi3lMET", &DPhi3lMET, &b_DPhi3lMET);
   fChain->SetBranchAddress("DPhi3lMET_up", &DPhi3lMET_up, &b_DPhi3lMET_up);
   fChain->SetBranchAddress("DPhi3lMET_dn", &DPhi3lMET_dn, &b_DPhi3lMET_dn);
   fChain->SetBranchAddress("DPhi3lMET_jer", &DPhi3lMET_jer, &b_DPhi3lMET_jer);
   fChain->SetBranchAddress("DPhi3lMET_jerup", &DPhi3lMET_jerup, &b_DPhi3lMET_jerup);
   fChain->SetBranchAddress("DPhi3lMET_jerdn", &DPhi3lMET_jerdn, &b_DPhi3lMET_jerdn);
   fChain->SetBranchAddress("DPhi3lMET_gen", &DPhi3lMET_gen, &b_DPhi3lMET_gen);
   fChain->SetBranchAddress("MTmax", &MTmax, &b_MTmax);
   fChain->SetBranchAddress("MTmax_up", &MTmax_up, &b_MTmax_up);
   fChain->SetBranchAddress("MTmax_dn", &MTmax_dn, &b_MTmax_dn);
   fChain->SetBranchAddress("MTmax_jer", &MTmax_jer, &b_MTmax_jer);
   fChain->SetBranchAddress("MTmax_jerup", &MTmax_jerup, &b_MTmax_jerup);
   fChain->SetBranchAddress("MTmax_jerdn", &MTmax_jerdn, &b_MTmax_jerdn);
   fChain->SetBranchAddress("MTmax_gen", &MTmax_gen, &b_MTmax_gen);
   fChain->SetBranchAddress("MTmin", &MTmin, &b_MTmin);
   fChain->SetBranchAddress("MTmin_up", &MTmin_up, &b_MTmin_up);
   fChain->SetBranchAddress("MTmin_dn", &MTmin_dn, &b_MTmin_dn);
   fChain->SetBranchAddress("MTmin_jer", &MTmin_jer, &b_MTmin_jer);
   fChain->SetBranchAddress("MTmin_jerup", &MTmin_jerup, &b_MTmin_jerup);
   fChain->SetBranchAddress("MTmin_jerdn", &MTmin_jerdn, &b_MTmin_jerdn);
   fChain->SetBranchAddress("MTmin_gen", &MTmin_gen, &b_MTmin_gen);
   fChain->SetBranchAddress("MT3rd", &MT3rd, &b_MT3rd);
   fChain->SetBranchAddress("MT3rd_up", &MT3rd_up, &b_MT3rd_up);
   fChain->SetBranchAddress("MT3rd_dn", &MT3rd_dn, &b_MT3rd_dn);
   fChain->SetBranchAddress("MT3rd_jer", &MT3rd_jer, &b_MT3rd_jer);
   fChain->SetBranchAddress("MT3rd_jerup", &MT3rd_jerup, &b_MT3rd_jerup);
   fChain->SetBranchAddress("MT3rd_jerdn", &MT3rd_jerdn, &b_MT3rd_jerdn);
   fChain->SetBranchAddress("MT3rd_gen", &MT3rd_gen, &b_MT3rd_gen);
   fChain->SetBranchAddress("MTmax3L", &MTmax3L, &b_MTmax3L);
   fChain->SetBranchAddress("MTmax3L_up", &MTmax3L_up, &b_MTmax3L_up);
   fChain->SetBranchAddress("MTmax3L_dn", &MTmax3L_dn, &b_MTmax3L_dn);
   fChain->SetBranchAddress("MTmax3L_jer", &MTmax3L_jer, &b_MTmax3L_jer);
   fChain->SetBranchAddress("MTmax3L_jerup", &MTmax3L_jerup, &b_MTmax3L_jerup);
   fChain->SetBranchAddress("MTmax3L_jerdn", &MTmax3L_jerdn, &b_MTmax3L_jerdn);
   fChain->SetBranchAddress("MTmax3L_gen", &MTmax3L_gen, &b_MTmax3L_gen);
   fChain->SetBranchAddress("passSSee", &passSSee, &b_passSSee);
   fChain->SetBranchAddress("passSSem", &passSSem, &b_passSSem);
   fChain->SetBranchAddress("passSSmm", &passSSmm, &b_passSSmm);
   fChain->SetBranchAddress("lep_idx0_SS", &lep_idx0_SS, &b_lep_idx0_SS);
   fChain->SetBranchAddress("lep_idx1_SS", &lep_idx1_SS, &b_lep_idx1_SS);
   fChain->SetBranchAddress("bkgtype", &bkgtype, &b_bkgtype);
   fChain->SetBranchAddress("vetophoton", &vetophoton, &b_vetophoton);
   fChain->SetBranchAddress("purewgt", &purewgt, &b_purewgt);
   fChain->SetBranchAddress("purewgt_up", &purewgt_up, &b_purewgt_up);
   fChain->SetBranchAddress("purewgt_dn", &purewgt_dn, &b_purewgt_dn);
   fChain->SetBranchAddress("ffwgt", &ffwgt, &b_ffwgt);
   fChain->SetBranchAddress("ffwgt_up", &ffwgt_up, &b_ffwgt_up);
   fChain->SetBranchAddress("ffwgt_dn", &ffwgt_dn, &b_ffwgt_dn);
   fChain->SetBranchAddress("ffwgt_el_up", &ffwgt_el_up, &b_ffwgt_el_up);
   fChain->SetBranchAddress("ffwgt_el_dn", &ffwgt_el_dn, &b_ffwgt_el_dn);
   fChain->SetBranchAddress("ffwgt_mu_up", &ffwgt_mu_up, &b_ffwgt_mu_up);
   fChain->SetBranchAddress("ffwgt_mu_dn", &ffwgt_mu_dn, &b_ffwgt_mu_dn);
   fChain->SetBranchAddress("ffwgt_closure_up", &ffwgt_closure_up, &b_ffwgt_closure_up);
   fChain->SetBranchAddress("ffwgt_closure_dn", &ffwgt_closure_dn, &b_ffwgt_closure_dn);
   fChain->SetBranchAddress("ffwgt_closure_el_up", &ffwgt_closure_el_up, &b_ffwgt_closure_el_up);
   fChain->SetBranchAddress("ffwgt_closure_el_dn", &ffwgt_closure_el_dn, &b_ffwgt_closure_el_dn);
   fChain->SetBranchAddress("ffwgt_closure_mu_up", &ffwgt_closure_mu_up, &b_ffwgt_closure_mu_up);
   fChain->SetBranchAddress("ffwgt_closure_mu_dn", &ffwgt_closure_mu_dn, &b_ffwgt_closure_mu_dn);
   fChain->SetBranchAddress("ffwgt_full_up", &ffwgt_full_up, &b_ffwgt_full_up);
   fChain->SetBranchAddress("ffwgt_full_dn", &ffwgt_full_dn, &b_ffwgt_full_dn);
   fChain->SetBranchAddress("ffwgtqcd", &ffwgtqcd, &b_ffwgtqcd);
   fChain->SetBranchAddress("ffwgtqcd_up", &ffwgtqcd_up, &b_ffwgtqcd_up);
   fChain->SetBranchAddress("ffwgtqcd_dn", &ffwgtqcd_dn, &b_ffwgtqcd_dn);
   fChain->SetBranchAddress("lepsf", &lepsf, &b_lepsf);
   fChain->SetBranchAddress("lepsf_up", &lepsf_up, &b_lepsf_up);
   fChain->SetBranchAddress("lepsf_dn", &lepsf_dn, &b_lepsf_dn);
   fChain->SetBranchAddress("trigeff", &trigeff, &b_trigeff);
   fChain->SetBranchAddress("trigeff_up", &trigeff_up, &b_trigeff_up);
   fChain->SetBranchAddress("trigeff_dn", &trigeff_dn, &b_trigeff_dn);
   fChain->SetBranchAddress("trigsf", &trigsf, &b_trigsf);
   fChain->SetBranchAddress("trigsf_up", &trigsf_up, &b_trigsf_up);
   fChain->SetBranchAddress("trigsf_dn", &trigsf_dn, &b_trigsf_dn);
   fChain->SetBranchAddress("musmear_sf", &musmear_sf, &b_musmear_sf);
   fChain->SetBranchAddress("higgsdecay", &higgsdecay, &b_higgsdecay);
   fChain->SetBranchAddress("nlep", &nlep, &b_nlep);
   fChain->SetBranchAddress("nquark", &nquark, &b_nquark);
   fChain->SetBranchAddress("wa_id", &wa_id, &b_wa_id);
//    fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_higgs_p4_fCoordinates_fX);
//    fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_higgs_p4_fCoordinates_fY);
//    fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_higgs_p4_fCoordinates_fZ);
//    fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_higgs_p4_fCoordinates_fT);
   fChain->SetBranchAddress("decay_p4", &decay_p4_, &b_decay_p4_);
   fChain->SetBranchAddress("decay_p4.fCoordinates.fX", decay_p4_fCoordinates_fX, &b_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("decay_p4.fCoordinates.fY", decay_p4_fCoordinates_fY, &b_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("decay_p4.fCoordinates.fZ", decay_p4_fCoordinates_fZ, &b_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("decay_p4.fCoordinates.fT", decay_p4_fCoordinates_fT, &b_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("decay_id", &decay_id, &b_decay_id);
   fChain->SetBranchAddress("decay_isstar", &decay_isstar, &b_decay_isstar);
   fChain->SetBranchAddress("lepton_p4", &lepton_p4_, &b_lepton_p4_);
   fChain->SetBranchAddress("lepton_p4.fCoordinates.fX", lepton_p4_fCoordinates_fX, &b_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("lepton_p4.fCoordinates.fY", lepton_p4_fCoordinates_fY, &b_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("lepton_p4.fCoordinates.fZ", lepton_p4_fCoordinates_fZ, &b_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("lepton_p4.fCoordinates.fT", lepton_p4_fCoordinates_fT, &b_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("lepton_id", &lepton_id, &b_lepton_id);
   fChain->SetBranchAddress("lepton_isstar", &lepton_isstar, &b_lepton_isstar);
   fChain->SetBranchAddress("quark_p4", &quark_p4_, &b_quark_p4_);
   fChain->SetBranchAddress("quark_p4.fCoordinates.fX", quark_p4_fCoordinates_fX, &b_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("quark_p4.fCoordinates.fY", quark_p4_fCoordinates_fY, &b_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("quark_p4.fCoordinates.fZ", quark_p4_fCoordinates_fZ, &b_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("quark_p4.fCoordinates.fT", quark_p4_fCoordinates_fT, &b_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("quark_id", &quark_id, &b_quark_id);
   fChain->SetBranchAddress("quark_isstar", &quark_isstar, &b_quark_isstar);
   fChain->SetBranchAddress("lqq_max_dr", &lqq_max_dr, &b_lqq_max_dr);
   fChain->SetBranchAddress("lq0_dr", &lq0_dr, &b_lq0_dr);
   fChain->SetBranchAddress("lq1_dr", &lq1_dr, &b_lq1_dr);
   fChain->SetBranchAddress("qq_dr", &qq_dr, &b_qq_dr);
   fChain->SetBranchAddress("boosted0_decay_p4", &boosted0_decay_p4_, &b_boosted0_decay_p4_);
   fChain->SetBranchAddress("boosted0_decay_p4.fCoordinates.fX", boosted0_decay_p4_fCoordinates_fX, &b_boosted0_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted0_decay_p4.fCoordinates.fY", boosted0_decay_p4_fCoordinates_fY, &b_boosted0_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted0_decay_p4.fCoordinates.fZ", boosted0_decay_p4_fCoordinates_fZ, &b_boosted0_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted0_decay_p4.fCoordinates.fT", boosted0_decay_p4_fCoordinates_fT, &b_boosted0_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted0_decay_id", &boosted0_decay_id, &b_boosted0_decay_id);
   fChain->SetBranchAddress("boosted0_decay_isstar", &boosted0_decay_isstar, &b_boosted0_decay_isstar);
   fChain->SetBranchAddress("boosted0_decay_h_dr", &boosted0_decay_h_dr, &b_boosted0_decay_h_dr);
   fChain->SetBranchAddress("boosted0_decay_h_deta", &boosted0_decay_h_deta, &b_boosted0_decay_h_deta);
   fChain->SetBranchAddress("boosted0_decay_h_dphi", &boosted0_decay_h_dphi, &b_boosted0_decay_h_dphi);
   fChain->SetBranchAddress("boosted0_decay_h_deta_rotated", &boosted0_decay_h_deta_rotated, &b_boosted0_decay_h_deta_rotated);
   fChain->SetBranchAddress("boosted0_decay_h_dphi_rotated", &boosted0_decay_h_dphi_rotated, &b_boosted0_decay_h_dphi_rotated);
   fChain->SetBranchAddress("boosted0_lepton_p4", &boosted0_lepton_p4_, &b_boosted0_lepton_p4_);
   fChain->SetBranchAddress("boosted0_lepton_p4.fCoordinates.fX", boosted0_lepton_p4_fCoordinates_fX, &b_boosted0_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted0_lepton_p4.fCoordinates.fY", boosted0_lepton_p4_fCoordinates_fY, &b_boosted0_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted0_lepton_p4.fCoordinates.fZ", boosted0_lepton_p4_fCoordinates_fZ, &b_boosted0_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted0_lepton_p4.fCoordinates.fT", boosted0_lepton_p4_fCoordinates_fT, &b_boosted0_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted0_lepton_id", &boosted0_lepton_id, &b_boosted0_lepton_id);
   fChain->SetBranchAddress("boosted0_lepton_isstar", &boosted0_lepton_isstar, &b_boosted0_lepton_isstar);
   fChain->SetBranchAddress("boosted0_lepton_h_dr", &boosted0_lepton_h_dr, &b_boosted0_lepton_h_dr);
   fChain->SetBranchAddress("boosted0_lepton_h_deta", &boosted0_lepton_h_deta, &b_boosted0_lepton_h_deta);
   fChain->SetBranchAddress("boosted0_lepton_h_dphi", &boosted0_lepton_h_dphi, &b_boosted0_lepton_h_dphi);
   fChain->SetBranchAddress("boosted0_lepton_h_deta_rotated", &boosted0_lepton_h_deta_rotated, &b_boosted0_lepton_h_deta_rotated);
   fChain->SetBranchAddress("boosted0_lepton_h_dphi_rotated", &boosted0_lepton_h_dphi_rotated, &b_boosted0_lepton_h_dphi_rotated);
   fChain->SetBranchAddress("boosted0_quark_p4", &boosted0_quark_p4_, &b_boosted0_quark_p4_);
   fChain->SetBranchAddress("boosted0_quark_p4.fCoordinates.fX", boosted0_quark_p4_fCoordinates_fX, &b_boosted0_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted0_quark_p4.fCoordinates.fY", boosted0_quark_p4_fCoordinates_fY, &b_boosted0_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted0_quark_p4.fCoordinates.fZ", boosted0_quark_p4_fCoordinates_fZ, &b_boosted0_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted0_quark_p4.fCoordinates.fT", boosted0_quark_p4_fCoordinates_fT, &b_boosted0_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted0_quark_id", &boosted0_quark_id, &b_boosted0_quark_id);
   fChain->SetBranchAddress("boosted0_quark_isstar", &boosted0_quark_isstar, &b_boosted0_quark_isstar);
   fChain->SetBranchAddress("boosted0_quark_h_dr", &boosted0_quark_h_dr, &b_boosted0_quark_h_dr);
   fChain->SetBranchAddress("boosted0_quark_h_deta", &boosted0_quark_h_deta, &b_boosted0_quark_h_deta);
   fChain->SetBranchAddress("boosted0_quark_h_dphi", &boosted0_quark_h_dphi, &b_boosted0_quark_h_dphi);
   fChain->SetBranchAddress("boosted0_quark_h_deta_rotated", &boosted0_quark_h_deta_rotated, &b_boosted0_quark_h_deta_rotated);
   fChain->SetBranchAddress("boosted0_quark_h_dphi_rotated", &boosted0_quark_h_dphi_rotated, &b_boosted0_quark_h_dphi_rotated);
   fChain->SetBranchAddress("boosted250_decay_p4", &boosted250_decay_p4_, &b_boosted250_decay_p4_);
   fChain->SetBranchAddress("boosted250_decay_p4.fCoordinates.fX", boosted250_decay_p4_fCoordinates_fX, &b_boosted250_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted250_decay_p4.fCoordinates.fY", boosted250_decay_p4_fCoordinates_fY, &b_boosted250_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted250_decay_p4.fCoordinates.fZ", boosted250_decay_p4_fCoordinates_fZ, &b_boosted250_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted250_decay_p4.fCoordinates.fT", boosted250_decay_p4_fCoordinates_fT, &b_boosted250_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted250_decay_id", &boosted250_decay_id, &b_boosted250_decay_id);
   fChain->SetBranchAddress("boosted250_decay_isstar", &boosted250_decay_isstar, &b_boosted250_decay_isstar);
   fChain->SetBranchAddress("boosted250_decay_h_dr", &boosted250_decay_h_dr, &b_boosted250_decay_h_dr);
   fChain->SetBranchAddress("boosted250_decay_h_deta", &boosted250_decay_h_deta, &b_boosted250_decay_h_deta);
   fChain->SetBranchAddress("boosted250_decay_h_dphi", &boosted250_decay_h_dphi, &b_boosted250_decay_h_dphi);
   fChain->SetBranchAddress("boosted250_decay_h_deta_rotated", &boosted250_decay_h_deta_rotated, &b_boosted250_decay_h_deta_rotated);
   fChain->SetBranchAddress("boosted250_decay_h_dphi_rotated", &boosted250_decay_h_dphi_rotated, &b_boosted250_decay_h_dphi_rotated);
   fChain->SetBranchAddress("boosted250_lepton_p4", &boosted250_lepton_p4_, &b_boosted250_lepton_p4_);
   fChain->SetBranchAddress("boosted250_lepton_p4.fCoordinates.fX", boosted250_lepton_p4_fCoordinates_fX, &b_boosted250_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted250_lepton_p4.fCoordinates.fY", boosted250_lepton_p4_fCoordinates_fY, &b_boosted250_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted250_lepton_p4.fCoordinates.fZ", boosted250_lepton_p4_fCoordinates_fZ, &b_boosted250_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted250_lepton_p4.fCoordinates.fT", boosted250_lepton_p4_fCoordinates_fT, &b_boosted250_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted250_lepton_id", &boosted250_lepton_id, &b_boosted250_lepton_id);
   fChain->SetBranchAddress("boosted250_lepton_isstar", &boosted250_lepton_isstar, &b_boosted250_lepton_isstar);
   fChain->SetBranchAddress("boosted250_lepton_h_dr", &boosted250_lepton_h_dr, &b_boosted250_lepton_h_dr);
   fChain->SetBranchAddress("boosted250_lepton_h_deta", &boosted250_lepton_h_deta, &b_boosted250_lepton_h_deta);
   fChain->SetBranchAddress("boosted250_lepton_h_dphi", &boosted250_lepton_h_dphi, &b_boosted250_lepton_h_dphi);
   fChain->SetBranchAddress("boosted250_lepton_h_deta_rotated", &boosted250_lepton_h_deta_rotated, &b_boosted250_lepton_h_deta_rotated);
   fChain->SetBranchAddress("boosted250_lepton_h_dphi_rotated", &boosted250_lepton_h_dphi_rotated, &b_boosted250_lepton_h_dphi_rotated);
   fChain->SetBranchAddress("boosted250_quark_p4", &boosted250_quark_p4_, &b_boosted250_quark_p4_);
   fChain->SetBranchAddress("boosted250_quark_p4.fCoordinates.fX", boosted250_quark_p4_fCoordinates_fX, &b_boosted250_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted250_quark_p4.fCoordinates.fY", boosted250_quark_p4_fCoordinates_fY, &b_boosted250_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted250_quark_p4.fCoordinates.fZ", boosted250_quark_p4_fCoordinates_fZ, &b_boosted250_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted250_quark_p4.fCoordinates.fT", boosted250_quark_p4_fCoordinates_fT, &b_boosted250_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted250_quark_id", &boosted250_quark_id, &b_boosted250_quark_id);
   fChain->SetBranchAddress("boosted250_quark_isstar", &boosted250_quark_isstar, &b_boosted250_quark_isstar);
   fChain->SetBranchAddress("boosted250_quark_h_dr", &boosted250_quark_h_dr, &b_boosted250_quark_h_dr);
   fChain->SetBranchAddress("boosted250_quark_h_deta", &boosted250_quark_h_deta, &b_boosted250_quark_h_deta);
   fChain->SetBranchAddress("boosted250_quark_h_dphi", &boosted250_quark_h_dphi, &b_boosted250_quark_h_dphi);
   fChain->SetBranchAddress("boosted250_quark_h_deta_rotated", &boosted250_quark_h_deta_rotated, &b_boosted250_quark_h_deta_rotated);
   fChain->SetBranchAddress("boosted250_quark_h_dphi_rotated", &boosted250_quark_h_dphi_rotated, &b_boosted250_quark_h_dphi_rotated);
   fChain->SetBranchAddress("boosted250_lqq_max_dr", &boosted250_lqq_max_dr, &b_boosted250_lqq_max_dr);
   fChain->SetBranchAddress("boosted250_lq0_dr", &boosted250_lq0_dr, &b_boosted250_lq0_dr);
   fChain->SetBranchAddress("boosted250_lq1_dr", &boosted250_lq1_dr, &b_boosted250_lq1_dr);
   fChain->SetBranchAddress("boosted250_qq_dr", &boosted250_qq_dr, &b_boosted250_qq_dr);
   fChain->SetBranchAddress("boosted500_decay_p4", &boosted500_decay_p4_, &b_boosted500_decay_p4_);
   fChain->SetBranchAddress("boosted500_decay_p4.fCoordinates.fX", boosted500_decay_p4_fCoordinates_fX, &b_boosted500_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted500_decay_p4.fCoordinates.fY", boosted500_decay_p4_fCoordinates_fY, &b_boosted500_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted500_decay_p4.fCoordinates.fZ", boosted500_decay_p4_fCoordinates_fZ, &b_boosted500_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted500_decay_p4.fCoordinates.fT", boosted500_decay_p4_fCoordinates_fT, &b_boosted500_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted500_decay_id", &boosted500_decay_id, &b_boosted500_decay_id);
   fChain->SetBranchAddress("boosted500_decay_isstar", &boosted500_decay_isstar, &b_boosted500_decay_isstar);
   fChain->SetBranchAddress("boosted500_decay_h_dr", &boosted500_decay_h_dr, &b_boosted500_decay_h_dr);
   fChain->SetBranchAddress("boosted500_decay_h_deta", &boosted500_decay_h_deta, &b_boosted500_decay_h_deta);
   fChain->SetBranchAddress("boosted500_decay_h_dphi", &boosted500_decay_h_dphi, &b_boosted500_decay_h_dphi);
   fChain->SetBranchAddress("boosted500_decay_h_deta_rotated", &boosted500_decay_h_deta_rotated, &b_boosted500_decay_h_deta_rotated);
   fChain->SetBranchAddress("boosted500_decay_h_dphi_rotated", &boosted500_decay_h_dphi_rotated, &b_boosted500_decay_h_dphi_rotated);
   fChain->SetBranchAddress("boosted500_lepton_p4", &boosted500_lepton_p4_, &b_boosted500_lepton_p4_);
   fChain->SetBranchAddress("boosted500_lepton_p4.fCoordinates.fX", boosted500_lepton_p4_fCoordinates_fX, &b_boosted500_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted500_lepton_p4.fCoordinates.fY", boosted500_lepton_p4_fCoordinates_fY, &b_boosted500_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted500_lepton_p4.fCoordinates.fZ", boosted500_lepton_p4_fCoordinates_fZ, &b_boosted500_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted500_lepton_p4.fCoordinates.fT", boosted500_lepton_p4_fCoordinates_fT, &b_boosted500_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted500_lepton_id", &boosted500_lepton_id, &b_boosted500_lepton_id);
   fChain->SetBranchAddress("boosted500_lepton_isstar", &boosted500_lepton_isstar, &b_boosted500_lepton_isstar);
   fChain->SetBranchAddress("boosted500_lepton_h_dr", &boosted500_lepton_h_dr, &b_boosted500_lepton_h_dr);
   fChain->SetBranchAddress("boosted500_lepton_h_deta", &boosted500_lepton_h_deta, &b_boosted500_lepton_h_deta);
   fChain->SetBranchAddress("boosted500_lepton_h_dphi", &boosted500_lepton_h_dphi, &b_boosted500_lepton_h_dphi);
   fChain->SetBranchAddress("boosted500_lepton_h_deta_rotated", &boosted500_lepton_h_deta_rotated, &b_boosted500_lepton_h_deta_rotated);
   fChain->SetBranchAddress("boosted500_lepton_h_dphi_rotated", &boosted500_lepton_h_dphi_rotated, &b_boosted500_lepton_h_dphi_rotated);
   fChain->SetBranchAddress("boosted500_quark_p4", &boosted500_quark_p4_, &b_boosted500_quark_p4_);
   fChain->SetBranchAddress("boosted500_quark_p4.fCoordinates.fX", boosted500_quark_p4_fCoordinates_fX, &b_boosted500_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted500_quark_p4.fCoordinates.fY", boosted500_quark_p4_fCoordinates_fY, &b_boosted500_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted500_quark_p4.fCoordinates.fZ", boosted500_quark_p4_fCoordinates_fZ, &b_boosted500_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted500_quark_p4.fCoordinates.fT", boosted500_quark_p4_fCoordinates_fT, &b_boosted500_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted500_quark_id", &boosted500_quark_id, &b_boosted500_quark_id);
   fChain->SetBranchAddress("boosted500_quark_isstar", &boosted500_quark_isstar, &b_boosted500_quark_isstar);
   fChain->SetBranchAddress("boosted500_quark_h_dr", &boosted500_quark_h_dr, &b_boosted500_quark_h_dr);
   fChain->SetBranchAddress("boosted500_quark_h_deta", &boosted500_quark_h_deta, &b_boosted500_quark_h_deta);
   fChain->SetBranchAddress("boosted500_quark_h_dphi", &boosted500_quark_h_dphi, &b_boosted500_quark_h_dphi);
   fChain->SetBranchAddress("boosted500_quark_h_deta_rotated", &boosted500_quark_h_deta_rotated, &b_boosted500_quark_h_deta_rotated);
   fChain->SetBranchAddress("boosted500_quark_h_dphi_rotated", &boosted500_quark_h_dphi_rotated, &b_boosted500_quark_h_dphi_rotated);
   fChain->SetBranchAddress("boosted500_lqq_max_dr", &boosted500_lqq_max_dr, &b_boosted500_lqq_max_dr);
   fChain->SetBranchAddress("boosted500_lq0_dr", &boosted500_lq0_dr, &b_boosted500_lq0_dr);
   fChain->SetBranchAddress("boosted500_lq1_dr", &boosted500_lq1_dr, &b_boosted500_lq1_dr);
   fChain->SetBranchAddress("boosted500_qq_dr", &boosted500_qq_dr, &b_boosted500_qq_dr);
   fChain->SetBranchAddress("boosted1000_decay_p4", &boosted1000_decay_p4_, &b_boosted1000_decay_p4_);
   fChain->SetBranchAddress("boosted1000_decay_p4.fCoordinates.fX", boosted1000_decay_p4_fCoordinates_fX, &b_boosted1000_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1000_decay_p4.fCoordinates.fY", boosted1000_decay_p4_fCoordinates_fY, &b_boosted1000_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1000_decay_p4.fCoordinates.fZ", boosted1000_decay_p4_fCoordinates_fZ, &b_boosted1000_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1000_decay_p4.fCoordinates.fT", boosted1000_decay_p4_fCoordinates_fT, &b_boosted1000_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1000_decay_id", &boosted1000_decay_id, &b_boosted1000_decay_id);
   fChain->SetBranchAddress("boosted1000_decay_isstar", &boosted1000_decay_isstar, &b_boosted1000_decay_isstar);
   fChain->SetBranchAddress("boosted1000_decay_h_dr", &boosted1000_decay_h_dr, &b_boosted1000_decay_h_dr);
   fChain->SetBranchAddress("boosted1000_decay_h_deta", &boosted1000_decay_h_deta, &b_boosted1000_decay_h_deta);
   fChain->SetBranchAddress("boosted1000_decay_h_dphi", &boosted1000_decay_h_dphi, &b_boosted1000_decay_h_dphi);
   fChain->SetBranchAddress("boosted1000_decay_h_deta_rotated", &boosted1000_decay_h_deta_rotated, &b_boosted1000_decay_h_deta_rotated);
   fChain->SetBranchAddress("boosted1000_decay_h_dphi_rotated", &boosted1000_decay_h_dphi_rotated, &b_boosted1000_decay_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1000_lepton_p4", &boosted1000_lepton_p4_, &b_boosted1000_lepton_p4_);
   fChain->SetBranchAddress("boosted1000_lepton_p4.fCoordinates.fX", boosted1000_lepton_p4_fCoordinates_fX, &b_boosted1000_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1000_lepton_p4.fCoordinates.fY", boosted1000_lepton_p4_fCoordinates_fY, &b_boosted1000_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1000_lepton_p4.fCoordinates.fZ", boosted1000_lepton_p4_fCoordinates_fZ, &b_boosted1000_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1000_lepton_p4.fCoordinates.fT", boosted1000_lepton_p4_fCoordinates_fT, &b_boosted1000_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1000_lepton_id", &boosted1000_lepton_id, &b_boosted1000_lepton_id);
   fChain->SetBranchAddress("boosted1000_lepton_isstar", &boosted1000_lepton_isstar, &b_boosted1000_lepton_isstar);
   fChain->SetBranchAddress("boosted1000_lepton_h_dr", &boosted1000_lepton_h_dr, &b_boosted1000_lepton_h_dr);
   fChain->SetBranchAddress("boosted1000_lepton_h_deta", &boosted1000_lepton_h_deta, &b_boosted1000_lepton_h_deta);
   fChain->SetBranchAddress("boosted1000_lepton_h_dphi", &boosted1000_lepton_h_dphi, &b_boosted1000_lepton_h_dphi);
   fChain->SetBranchAddress("boosted1000_lepton_h_deta_rotated", &boosted1000_lepton_h_deta_rotated, &b_boosted1000_lepton_h_deta_rotated);
   fChain->SetBranchAddress("boosted1000_lepton_h_dphi_rotated", &boosted1000_lepton_h_dphi_rotated, &b_boosted1000_lepton_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1000_quark_p4", &boosted1000_quark_p4_, &b_boosted1000_quark_p4_);
   fChain->SetBranchAddress("boosted1000_quark_p4.fCoordinates.fX", boosted1000_quark_p4_fCoordinates_fX, &b_boosted1000_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1000_quark_p4.fCoordinates.fY", boosted1000_quark_p4_fCoordinates_fY, &b_boosted1000_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1000_quark_p4.fCoordinates.fZ", boosted1000_quark_p4_fCoordinates_fZ, &b_boosted1000_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1000_quark_p4.fCoordinates.fT", boosted1000_quark_p4_fCoordinates_fT, &b_boosted1000_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1000_quark_id", &boosted1000_quark_id, &b_boosted1000_quark_id);
   fChain->SetBranchAddress("boosted1000_quark_isstar", &boosted1000_quark_isstar, &b_boosted1000_quark_isstar);
   fChain->SetBranchAddress("boosted1000_quark_h_dr", &boosted1000_quark_h_dr, &b_boosted1000_quark_h_dr);
   fChain->SetBranchAddress("boosted1000_quark_h_deta", &boosted1000_quark_h_deta, &b_boosted1000_quark_h_deta);
   fChain->SetBranchAddress("boosted1000_quark_h_dphi", &boosted1000_quark_h_dphi, &b_boosted1000_quark_h_dphi);
   fChain->SetBranchAddress("boosted1000_quark_h_deta_rotated", &boosted1000_quark_h_deta_rotated, &b_boosted1000_quark_h_deta_rotated);
   fChain->SetBranchAddress("boosted1000_quark_h_dphi_rotated", &boosted1000_quark_h_dphi_rotated, &b_boosted1000_quark_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1000_lqq_max_dr", &boosted1000_lqq_max_dr, &b_boosted1000_lqq_max_dr);
   fChain->SetBranchAddress("boosted1000_lq0_dr", &boosted1000_lq0_dr, &b_boosted1000_lq0_dr);
   fChain->SetBranchAddress("boosted1000_lq1_dr", &boosted1000_lq1_dr, &b_boosted1000_lq1_dr);
   fChain->SetBranchAddress("boosted1000_qq_dr", &boosted1000_qq_dr, &b_boosted1000_qq_dr);
   fChain->SetBranchAddress("boosted1500_decay_p4", &boosted1500_decay_p4_, &b_boosted1500_decay_p4_);
   fChain->SetBranchAddress("boosted1500_decay_p4.fCoordinates.fX", boosted1500_decay_p4_fCoordinates_fX, &b_boosted1500_decay_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1500_decay_p4.fCoordinates.fY", boosted1500_decay_p4_fCoordinates_fY, &b_boosted1500_decay_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1500_decay_p4.fCoordinates.fZ", boosted1500_decay_p4_fCoordinates_fZ, &b_boosted1500_decay_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1500_decay_p4.fCoordinates.fT", boosted1500_decay_p4_fCoordinates_fT, &b_boosted1500_decay_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1500_decay_id", &boosted1500_decay_id, &b_boosted1500_decay_id);
   fChain->SetBranchAddress("boosted1500_decay_isstar", &boosted1500_decay_isstar, &b_boosted1500_decay_isstar);
   fChain->SetBranchAddress("boosted1500_decay_h_dr", &boosted1500_decay_h_dr, &b_boosted1500_decay_h_dr);
   fChain->SetBranchAddress("boosted1500_decay_h_deta", &boosted1500_decay_h_deta, &b_boosted1500_decay_h_deta);
   fChain->SetBranchAddress("boosted1500_decay_h_dphi", &boosted1500_decay_h_dphi, &b_boosted1500_decay_h_dphi);
   fChain->SetBranchAddress("boosted1500_decay_h_deta_rotated", &boosted1500_decay_h_deta_rotated, &b_boosted1500_decay_h_deta_rotated);
   fChain->SetBranchAddress("boosted1500_decay_h_dphi_rotated", &boosted1500_decay_h_dphi_rotated, &b_boosted1500_decay_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1500_lepton_p4", &boosted1500_lepton_p4_, &b_boosted1500_lepton_p4_);
   fChain->SetBranchAddress("boosted1500_lepton_p4.fCoordinates.fX", boosted1500_lepton_p4_fCoordinates_fX, &b_boosted1500_lepton_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1500_lepton_p4.fCoordinates.fY", boosted1500_lepton_p4_fCoordinates_fY, &b_boosted1500_lepton_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1500_lepton_p4.fCoordinates.fZ", boosted1500_lepton_p4_fCoordinates_fZ, &b_boosted1500_lepton_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1500_lepton_p4.fCoordinates.fT", boosted1500_lepton_p4_fCoordinates_fT, &b_boosted1500_lepton_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1500_lepton_id", &boosted1500_lepton_id, &b_boosted1500_lepton_id);
   fChain->SetBranchAddress("boosted1500_lepton_isstar", &boosted1500_lepton_isstar, &b_boosted1500_lepton_isstar);
   fChain->SetBranchAddress("boosted1500_lepton_h_dr", &boosted1500_lepton_h_dr, &b_boosted1500_lepton_h_dr);
   fChain->SetBranchAddress("boosted1500_lepton_h_deta", &boosted1500_lepton_h_deta, &b_boosted1500_lepton_h_deta);
   fChain->SetBranchAddress("boosted1500_lepton_h_dphi", &boosted1500_lepton_h_dphi, &b_boosted1500_lepton_h_dphi);
   fChain->SetBranchAddress("boosted1500_lepton_h_deta_rotated", &boosted1500_lepton_h_deta_rotated, &b_boosted1500_lepton_h_deta_rotated);
   fChain->SetBranchAddress("boosted1500_lepton_h_dphi_rotated", &boosted1500_lepton_h_dphi_rotated, &b_boosted1500_lepton_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1500_quark_p4", &boosted1500_quark_p4_, &b_boosted1500_quark_p4_);
   fChain->SetBranchAddress("boosted1500_quark_p4.fCoordinates.fX", boosted1500_quark_p4_fCoordinates_fX, &b_boosted1500_quark_p4_fCoordinates_fX);
   fChain->SetBranchAddress("boosted1500_quark_p4.fCoordinates.fY", boosted1500_quark_p4_fCoordinates_fY, &b_boosted1500_quark_p4_fCoordinates_fY);
   fChain->SetBranchAddress("boosted1500_quark_p4.fCoordinates.fZ", boosted1500_quark_p4_fCoordinates_fZ, &b_boosted1500_quark_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("boosted1500_quark_p4.fCoordinates.fT", boosted1500_quark_p4_fCoordinates_fT, &b_boosted1500_quark_p4_fCoordinates_fT);
   fChain->SetBranchAddress("boosted1500_quark_id", &boosted1500_quark_id, &b_boosted1500_quark_id);
   fChain->SetBranchAddress("boosted1500_quark_isstar", &boosted1500_quark_isstar, &b_boosted1500_quark_isstar);
   fChain->SetBranchAddress("boosted1500_quark_h_dr", &boosted1500_quark_h_dr, &b_boosted1500_quark_h_dr);
   fChain->SetBranchAddress("boosted1500_quark_h_deta", &boosted1500_quark_h_deta, &b_boosted1500_quark_h_deta);
   fChain->SetBranchAddress("boosted1500_quark_h_dphi", &boosted1500_quark_h_dphi, &b_boosted1500_quark_h_dphi);
   fChain->SetBranchAddress("boosted1500_quark_h_deta_rotated", &boosted1500_quark_h_deta_rotated, &b_boosted1500_quark_h_deta_rotated);
   fChain->SetBranchAddress("boosted1500_quark_h_dphi_rotated", &boosted1500_quark_h_dphi_rotated, &b_boosted1500_quark_h_dphi_rotated);
   fChain->SetBranchAddress("boosted1500_lqq_max_dr", &boosted1500_lqq_max_dr, &b_boosted1500_lqq_max_dr);
   fChain->SetBranchAddress("boosted1500_lq0_dr", &boosted1500_lq0_dr, &b_boosted1500_lq0_dr);
   fChain->SetBranchAddress("boosted1500_lq1_dr", &boosted1500_lq1_dr, &b_boosted1500_lq1_dr);
   fChain->SetBranchAddress("boosted1500_qq_dr", &boosted1500_qq_dr, &b_boosted1500_qq_dr);
   fChain->SetBranchAddress("iswhwww", &iswhwww, &b_iswhwww);
   fChain->SetBranchAddress("www_channel", &www_channel, &b_www_channel);
   fChain->SetBranchAddress("has_tau", &has_tau, &b_has_tau);
   fChain->SetBranchAddress("w_p4", &w_p4_, &b_w_p4_);
   fChain->SetBranchAddress("w_p4.fCoordinates.fX", w_p4_fCoordinates_fX, &b_w_p4_fCoordinates_fX);
   fChain->SetBranchAddress("w_p4.fCoordinates.fY", w_p4_fCoordinates_fY, &b_w_p4_fCoordinates_fY);
   fChain->SetBranchAddress("w_p4.fCoordinates.fZ", w_p4_fCoordinates_fZ, &b_w_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("w_p4.fCoordinates.fT", w_p4_fCoordinates_fT, &b_w_p4_fCoordinates_fT);
   fChain->SetBranchAddress("w_islep", &w_islep, &b_w_islep);
   fChain->SetBranchAddress("w_isstar", &w_isstar, &b_w_isstar);
   fChain->SetBranchAddress("w_isH", &w_isH, &b_w_isH);
   fChain->SetBranchAddress("l_p4", &l_p4_, &b_l_p4_);
   fChain->SetBranchAddress("l_p4.fCoordinates.fX", l_p4_fCoordinates_fX, &b_l_p4_fCoordinates_fX);
   fChain->SetBranchAddress("l_p4.fCoordinates.fY", l_p4_fCoordinates_fY, &b_l_p4_fCoordinates_fY);
   fChain->SetBranchAddress("l_p4.fCoordinates.fZ", l_p4_fCoordinates_fZ, &b_l_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("l_p4.fCoordinates.fT", l_p4_fCoordinates_fT, &b_l_p4_fCoordinates_fT);
   fChain->SetBranchAddress("l_w_pt", &l_w_pt, &b_l_w_pt);
   fChain->SetBranchAddress("l_w_eta", &l_w_eta, &b_l_w_eta);
   fChain->SetBranchAddress("l_w_phi", &l_w_phi, &b_l_w_phi);
   fChain->SetBranchAddress("l_w_mass", &l_w_mass, &b_l_w_mass);
   fChain->SetBranchAddress("l_w_id", &l_w_id, &b_l_w_id);
   fChain->SetBranchAddress("l_isstar", &l_isstar, &b_l_isstar);
   fChain->SetBranchAddress("l_isH", &l_isH, &b_l_isH);
   fChain->SetBranchAddress("l_istau", &l_istau, &b_l_istau);
   fChain->SetBranchAddress("q_p4", &q_p4_, &b_q_p4_);
   fChain->SetBranchAddress("q_p4.fCoordinates.fX", q_p4_fCoordinates_fX, &b_q_p4_fCoordinates_fX);
   fChain->SetBranchAddress("q_p4.fCoordinates.fY", q_p4_fCoordinates_fY, &b_q_p4_fCoordinates_fY);
   fChain->SetBranchAddress("q_p4.fCoordinates.fZ", q_p4_fCoordinates_fZ, &b_q_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("q_p4.fCoordinates.fT", q_p4_fCoordinates_fT, &b_q_p4_fCoordinates_fT);
   fChain->SetBranchAddress("q_w_pt", &q_w_pt, &b_q_w_pt);
   fChain->SetBranchAddress("q_w_eta", &q_w_eta, &b_q_w_eta);
   fChain->SetBranchAddress("q_w_phi", &q_w_phi, &b_q_w_phi);
   fChain->SetBranchAddress("q_w_mass", &q_w_mass, &b_q_w_mass);
   fChain->SetBranchAddress("q_w_id", &q_w_id, &b_q_w_id);
   fChain->SetBranchAddress("q_isstar", &q_isstar, &b_q_isstar);
   fChain->SetBranchAddress("q_isH", &q_isH, &b_q_isH);
   fChain->SetBranchAddress("dRllSS", &dRllSS, &b_dRllSS);
   fChain->SetBranchAddress("dRqqSS", &dRqqSS, &b_dRqqSS);
   fChain->SetBranchAddress("DPhill_higgs", &DPhill_higgs, &b_DPhill_higgs);
   fChain->SetBranchAddress("Mll_higgs", &Mll_higgs, &b_Mll_higgs);
   fChain->SetBranchAddress("MT_higgs", &MT_higgs, &b_MT_higgs);
   Notify();
}

Bool_t t_www::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t_www::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t_www::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t_www_cxx
