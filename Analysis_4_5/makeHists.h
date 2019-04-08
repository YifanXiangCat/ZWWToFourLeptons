#ifndef _MAKE_HISTS_H_
#define _MAKE_HISTS_H_
#include <iostream>
class TH1D;
class TH2D;
class TH2F;
class TH3F;
class TFile;
class TProfile;
class TProfile2D;
class makeHists
{
 public:
  TFile *hf;
//*********Oiginal histograms***********//

    TH1D *plot_cut_tri_mll;
    TH1D *plot_cut_tri_deltaR;
    TH1D *plot_cut_tri_deltaPhi;

    TH1D *plot_cut_tri_n0_nsfos_mll;
    TH1D *plot_cut_tri_n0_nsfos_deltaR;
    TH1D *plot_cut_tri_n0_nsfos_deltaPhi;

    TH1D *plot_cut_tri_n0_nj_mll;
    TH1D *plot_cut_tri_n0_nj_deltaR;
    TH1D *plot_cut_tri_n0_nj_deltaPhi;

    TH1D *plot_cut_tri_n0_nb_mll;
    TH1D *plot_cut_tri_n0_nb_deltaR;
    TH1D *plot_cut_tri_n0_nb_deltaPhi;

    TH1D *plot_cut_tri_n0_dphi_mll;
    TH1D *plot_cut_tri_n0_dphi_deltaR;
    TH1D *plot_cut_tri_n0_dphi_deltaPhi;

    TH1D *plot_cut_tri_n0_metpt_mll;
    TH1D *plot_cut_tri_n0_metpt_deltaR;
    TH1D *plot_cut_tri_n0_metpt_deltaPhi;

    TH1D *plot_cut_tri_n0_mll3l_mll;
    TH1D *plot_cut_tri_n0_mll3l_deltaR;
    TH1D *plot_cut_tri_n0_mll3l_deltaPhi;

    TH1D *plot_cut_tri_n0_m3l_mll;
    TH1D *plot_cut_tri_n0_m3l_deltaR;
    TH1D *plot_cut_tri_n0_m3l_deltaPhi;

    TH1D *plot_cut_tri_n0_mee3l_mll;
    TH1D *plot_cut_tri_n0_mee3l_deltaR;
    TH1D *plot_cut_tri_n0_mee3l_deltaPhi;

    TH1D *plot_cut_tri_n0_mtmax_mll;
    TH1D *plot_cut_tri_n0_mtmax_deltaR;
    TH1D *plot_cut_tri_n0_mtmax_deltaPhi;



    TH1D *plot_cut_tri_n1_nsfos_mll;
    TH1D *plot_cut_tri_n1_nsfos_deltaR;
    TH1D *plot_cut_tri_n1_nsfos_deltaPhi;

    TH1D *plot_cut_tri_n1_nj_mll;
    TH1D *plot_cut_tri_n1_nj_deltaR;
    TH1D *plot_cut_tri_n1_nj_deltaPhi;

    TH1D *plot_cut_tri_n1_nb_mll;
    TH1D *plot_cut_tri_n1_nb_deltaR;
    TH1D *plot_cut_tri_n1_nb_deltaPhi;

    TH1D *plot_cut_tri_n1_pt3l_mll;
    TH1D *plot_cut_tri_n1_pt3l_deltaR;
    TH1D *plot_cut_tri_n1_pt3l_deltaPhi;

    TH1D *plot_cut_tri_n1_dphi_mll;
    TH1D *plot_cut_tri_n1_dphi_deltaR;
    TH1D *plot_cut_tri_n1_dphi_deltaPhi;

    TH1D *plot_cut_tri_n1_metpt_mll;
    TH1D *plot_cut_tri_n1_metpt_deltaR;
    TH1D *plot_cut_tri_n1_metpt_deltaPhi;

    TH1D *plot_cut_tri_n1_mll3l_mll;
    TH1D *plot_cut_tri_n1_mll3l_deltaR;
    TH1D *plot_cut_tri_n1_mll3l_deltaPhi;

    TH1D *plot_cut_tri_n1_m3l_mll;
    TH1D *plot_cut_tri_n1_m3l_deltaR;
    TH1D *plot_cut_tri_n1_m3l_deltaPhi;

    TH1D *plot_cut_tri_n1_nsfosinz_mll;
    TH1D *plot_cut_tri_n1_nsfosinz_deltaR;
    TH1D *plot_cut_tri_n1_nsfosinz_deltaPhi;

    TH1D *plot_cut_tri_n1_mt3rd_mll;
    TH1D *plot_cut_tri_n1_mt3rd_deltaR;
    TH1D *plot_cut_tri_n1_mt3rd_deltaPhi;



    TH1D *plot_cut_tri_n2_nsfos_mll;
    TH1D *plot_cut_tri_n2_nsfos_deltaR;
    TH1D *plot_cut_tri_n2_nsfos_deltaPhi;

    TH1D *plot_cut_tri_n2_nj_mll;
    TH1D *plot_cut_tri_n2_nj_deltaR;
    TH1D *plot_cut_tri_n2_nj_deltaPhi;

    TH1D *plot_cut_tri_n2_nb_mll;
    TH1D *plot_cut_tri_n2_nb_deltaR;
    TH1D *plot_cut_tri_n2_nb_deltaPhi;

    TH1D *plot_cut_tri_n2_pt3l_mll;
    TH1D *plot_cut_tri_n2_pt3l_deltaR;
    TH1D *plot_cut_tri_n2_pt3l_deltaPhi;

    TH1D *plot_cut_tri_n2_dphi_mll;
    TH1D *plot_cut_tri_n2_dphi_deltaR;
    TH1D *plot_cut_tri_n2_dphi_deltaPhi;

    TH1D *plot_cut_tri_n2_metpt_mll;
    TH1D *plot_cut_tri_n2_metpt_deltaR;
    TH1D *plot_cut_tri_n2_metpt_deltaPhi;

    TH1D *plot_cut_tri_n2_mll3l_mll;
    TH1D *plot_cut_tri_n2_mll3l_deltaR;
    TH1D *plot_cut_tri_n2_mll3l_deltaPhi;

    TH1D *plot_cut_tri_n2_m3l_mll;
    TH1D *plot_cut_tri_n2_m3l_deltaR;
    TH1D *plot_cut_tri_n2_m3l_deltaPhi;

    TH1D *plot_cut_tri_n2_nsfosinz_mll;
    TH1D *plot_cut_tri_n2_nsfosinz_deltaR;
    TH1D *plot_cut_tri_n2_nsfosinz_deltaPhi;





//  TH1D *plot_zboson_mass_seed_muon[16];
//  TH1D *plot_zboson_mass_seed_antimuon[16];
//  TH1D *plot_zboson_mass_forward_muon[16];
//  TH1D *plot_zboson_mass_forward_antimuon[16];

//  TH1D *plot_zboson_mass;
//  TH1D *plot_zboson_mass_forward;
//  TH1D *plot_zboson_mass_backward;

//*************************************//
  void bookHists(const char* fName);
  void saveHists();

};

#endif
