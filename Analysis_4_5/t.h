//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  5 15:54:31 2019 by ROOT version 6.14/04
// from TTree t/All events
// found on file: wwz_amcatnlo_1.root
//////////////////////////////////////////////////////////

#ifndef t_h
#define t_h

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

class t {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxlep_p4 = 4;
   static constexpr Int_t kMaxjets_p4 = 15;

   // Declaration of leaf types
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
   Int_t           lep_p4_;
   Float_t         lep_p4_fCoordinates_fX[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fY[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fZ[kMaxlep_p4];   //[lep_p4_]
   Float_t         lep_p4_fCoordinates_fT[kMaxlep_p4];   //[lep_p4_]
   vector<float>   *lep_pt;
   vector<float>   *lep_eta;
   vector<float>   *lep_phi;
   vector<float>   *lep_energy;
   vector<float>   *lep_mva;
   vector<float>   *lep_relIso03EA;
   vector<float>   *lep_relIso03EAwLep;
   vector<float>   *lep_ip3d;
   vector<float>   *lep_sip3d;
   vector<float>   *lep_dxy;
   vector<float>   *lep_dz;
   vector<int>     *lep_mc_id;
   vector<int>     *lep_motherIdv2;
   vector<int>     *lep_idx;
   vector<int>     *lep_id;
   vector<int>     *lep_isTightPOG;
   vector<int>     *lep_isMediumPOG;
 //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *met_p4;
   Float_t         fCoordinates_fX;
   Float_t         fCoordinates_fY;
   Float_t         fCoordinates_fZ;
   Float_t         fCoordinates_fT;
   Float_t         met_pt;
   Float_t         met_phi;
   Float_t         met_up_pt;
   Float_t         met_up_phi;
   Float_t         met_dn_pt;
   Float_t         met_dn_phi;
   Float_t         met_gen_pt;
   Float_t         met_gen_phi;
   Int_t           jets_p4_;
   Float_t         jets_p4_fCoordinates_fX[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fY[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fZ[kMaxjets_p4];   //[jets_p4_]
   Float_t         jets_p4_fCoordinates_fT[kMaxjets_p4];   //[jets_p4_]
   vector<float>   *jets_pt;
   vector<float>   *jets_eta;
   vector<float>   *jets_phi;
   vector<float>   *jets_mass;
   Int_t           nj;
   Int_t           nb;
   Int_t           nbmed;
   Float_t         ht;
   Float_t         weight_btagsf;
   Float_t         weight_btagsf_heavy_DN;
   Float_t         weight_btagsf_heavy_UP;
   Float_t         weight_btagsf_light_DN;
   Float_t         weight_btagsf_light_UP;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_evt_scale1fb;   //!
   TBranch        *b_genps_weight;   //!
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
   TBranch        *b_HLT_DoubleMu;   //!
   TBranch        *b_HLT_DoubleEl;   //!
   TBranch        *b_HLT_MuEG;   //!
   TBranch        *b_pass_duplicate_ee_em_mm;   //!
   TBranch        *b_pass_duplicate_mm_em_ee;   //!
   TBranch        *b_gen_ht;   //!
   TBranch        *b_firstgoodvertex;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_nTrueInt;   //!
   TBranch        *b_lep_p4_;   //!
   TBranch        *b_lep_p4_fCoordinates_fX;   //!
   TBranch        *b_lep_p4_fCoordinates_fY;   //!
   TBranch        *b_lep_p4_fCoordinates_fZ;   //!
   TBranch        *b_lep_p4_fCoordinates_fT;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_energy;   //!
   TBranch        *b_lep_mva;   //!
   TBranch        *b_lep_relIso03EA;   //!
   TBranch        *b_lep_relIso03EAwLep;   //!
   TBranch        *b_lep_ip3d;   //!
   TBranch        *b_lep_sip3d;   //!
   TBranch        *b_lep_dxy;   //!
   TBranch        *b_lep_dz;   //!
   TBranch        *b_lep_mc_id;   //!
   TBranch        *b_lep_motherIdv2;   //!
   TBranch        *b_lep_idx;   //!
   TBranch        *b_lep_id;   //!
   TBranch        *b_lep_isTightPOG;   //!
   TBranch        *b_lep_isMediumPOG;   //!
   TBranch        *b_met_p4_fCoordinates_fX;   //!
   TBranch        *b_met_p4_fCoordinates_fY;   //!
   TBranch        *b_met_p4_fCoordinates_fZ;   //!
   TBranch        *b_met_p4_fCoordinates_fT;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_up_pt;   //!
   TBranch        *b_met_up_phi;   //!
   TBranch        *b_met_dn_pt;   //!
   TBranch        *b_met_dn_phi;   //!
   TBranch        *b_met_gen_pt;   //!
   TBranch        *b_met_gen_phi;   //!
   TBranch        *b_jets_p4_;   //!
   TBranch        *b_jets_p4_fCoordinates_fX;   //!
   TBranch        *b_jets_p4_fCoordinates_fY;   //!
   TBranch        *b_jets_p4_fCoordinates_fZ;   //!
   TBranch        *b_jets_p4_fCoordinates_fT;   //!
   TBranch        *b_jets_pt;   //!
   TBranch        *b_jets_eta;   //!
   TBranch        *b_jets_phi;   //!
   TBranch        *b_jets_mass;   //!
   TBranch        *b_nj;   //!
   TBranch        *b_nb;   //!
   TBranch        *b_nbmed;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_weight_btagsf;   //!
   TBranch        *b_weight_btagsf_heavy_DN;   //!
   TBranch        *b_weight_btagsf_heavy_UP;   //!
   TBranch        *b_weight_btagsf_light_DN;   //!
   TBranch        *b_weight_btagsf_light_UP;   //!

   t(TTree *tree=0);
   virtual ~t();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef t_cxx
t::t(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("wwz_amcatnlo_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("wwz_amcatnlo_1.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

t::~t()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t t::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t t::LoadTree(Long64_t entry)
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

void t::Init(TTree *tree)
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
   lep_energy = 0;
   lep_mva = 0;
   lep_relIso03EA = 0;
   lep_relIso03EAwLep = 0;
   lep_ip3d = 0;
   lep_sip3d = 0;
   lep_dxy = 0;
   lep_dz = 0;
   lep_mc_id = 0;
   lep_motherIdv2 = 0;
   lep_idx = 0;
   lep_id = 0;
   lep_isTightPOG = 0;
   lep_isMediumPOG = 0;
   jets_pt = 0;
   jets_eta = 0;
   jets_phi = 0;
   jets_mass = 0;
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
   fChain->SetBranchAddress("genps_weight", &genps_weight, &b_genps_weight);
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
   fChain->SetBranchAddress("HLT_DoubleMu", &HLT_DoubleMu, &b_HLT_DoubleMu);
   fChain->SetBranchAddress("HLT_DoubleEl", &HLT_DoubleEl, &b_HLT_DoubleEl);
   fChain->SetBranchAddress("HLT_MuEG", &HLT_MuEG, &b_HLT_MuEG);
   fChain->SetBranchAddress("pass_duplicate_ee_em_mm", &pass_duplicate_ee_em_mm, &b_pass_duplicate_ee_em_mm);
   fChain->SetBranchAddress("pass_duplicate_mm_em_ee", &pass_duplicate_mm_em_ee, &b_pass_duplicate_mm_em_ee);
   fChain->SetBranchAddress("gen_ht", &gen_ht, &b_gen_ht);
   fChain->SetBranchAddress("firstgoodvertex", &firstgoodvertex, &b_firstgoodvertex);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("nTrueInt", &nTrueInt, &b_nTrueInt);
   fChain->SetBranchAddress("lep_p4", &lep_p4_, &b_lep_p4_);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fX", lep_p4_fCoordinates_fX, &b_lep_p4_fCoordinates_fX);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fY", lep_p4_fCoordinates_fY, &b_lep_p4_fCoordinates_fY);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fZ", lep_p4_fCoordinates_fZ, &b_lep_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("lep_p4.fCoordinates.fT", lep_p4_fCoordinates_fT, &b_lep_p4_fCoordinates_fT);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_energy", &lep_energy, &b_lep_energy);
   fChain->SetBranchAddress("lep_mva", &lep_mva, &b_lep_mva);
   fChain->SetBranchAddress("lep_relIso03EA", &lep_relIso03EA, &b_lep_relIso03EA);
   fChain->SetBranchAddress("lep_relIso03EAwLep", &lep_relIso03EAwLep, &b_lep_relIso03EAwLep);
   fChain->SetBranchAddress("lep_ip3d", &lep_ip3d, &b_lep_ip3d);
   fChain->SetBranchAddress("lep_sip3d", &lep_sip3d, &b_lep_sip3d);
   fChain->SetBranchAddress("lep_dxy", &lep_dxy, &b_lep_dxy);
   fChain->SetBranchAddress("lep_dz", &lep_dz, &b_lep_dz);
   fChain->SetBranchAddress("lep_mc_id", &lep_mc_id, &b_lep_mc_id);
   fChain->SetBranchAddress("lep_motherIdv2", &lep_motherIdv2, &b_lep_motherIdv2);
   fChain->SetBranchAddress("lep_idx", &lep_idx, &b_lep_idx);
   fChain->SetBranchAddress("lep_id", &lep_id, &b_lep_id);
   fChain->SetBranchAddress("lep_isTightPOG", &lep_isTightPOG, &b_lep_isTightPOG);
   fChain->SetBranchAddress("lep_isMediumPOG", &lep_isMediumPOG, &b_lep_isMediumPOG);
   fChain->SetBranchAddress("fCoordinates.fX", &fCoordinates_fX, &b_met_p4_fCoordinates_fX);
   fChain->SetBranchAddress("fCoordinates.fY", &fCoordinates_fY, &b_met_p4_fCoordinates_fY);
   fChain->SetBranchAddress("fCoordinates.fZ", &fCoordinates_fZ, &b_met_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("fCoordinates.fT", &fCoordinates_fT, &b_met_p4_fCoordinates_fT);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_up_pt", &met_up_pt, &b_met_up_pt);
   fChain->SetBranchAddress("met_up_phi", &met_up_phi, &b_met_up_phi);
   fChain->SetBranchAddress("met_dn_pt", &met_dn_pt, &b_met_dn_pt);
   fChain->SetBranchAddress("met_dn_phi", &met_dn_phi, &b_met_dn_phi);
   fChain->SetBranchAddress("met_gen_pt", &met_gen_pt, &b_met_gen_pt);
   fChain->SetBranchAddress("met_gen_phi", &met_gen_phi, &b_met_gen_phi);
   fChain->SetBranchAddress("jets_p4", &jets_p4_, &b_jets_p4_);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fX", jets_p4_fCoordinates_fX, &b_jets_p4_fCoordinates_fX);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fY", jets_p4_fCoordinates_fY, &b_jets_p4_fCoordinates_fY);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fZ", jets_p4_fCoordinates_fZ, &b_jets_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("jets_p4.fCoordinates.fT", jets_p4_fCoordinates_fT, &b_jets_p4_fCoordinates_fT);
   fChain->SetBranchAddress("jets_pt", &jets_pt, &b_jets_pt);
   fChain->SetBranchAddress("jets_eta", &jets_eta, &b_jets_eta);
   fChain->SetBranchAddress("jets_phi", &jets_phi, &b_jets_phi);
   fChain->SetBranchAddress("jets_mass", &jets_mass, &b_jets_mass);
   fChain->SetBranchAddress("nj", &nj, &b_nj);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("nbmed", &nbmed, &b_nbmed);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("weight_btagsf", &weight_btagsf, &b_weight_btagsf);
   fChain->SetBranchAddress("weight_btagsf_heavy_DN", &weight_btagsf_heavy_DN, &b_weight_btagsf_heavy_DN);
   fChain->SetBranchAddress("weight_btagsf_heavy_UP", &weight_btagsf_heavy_UP, &b_weight_btagsf_heavy_UP);
   fChain->SetBranchAddress("weight_btagsf_light_DN", &weight_btagsf_light_DN, &b_weight_btagsf_light_DN);
   fChain->SetBranchAddress("weight_btagsf_light_UP", &weight_btagsf_light_UP, &b_weight_btagsf_light_UP);
   Notify();
}

Bool_t t::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void t::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t t::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef t_cxx
