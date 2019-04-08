
#include "makeHists.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
  
using namespace std;
void makeHists::bookHists(const char* fName)
{
	hf = new TFile(fName,"RECREATE");
        char nameHLT[100];


	double Rmllh=100;
	double Rmlll=0;
	int    Bmll=180;
	double RdeltaRh=4;
	double RdeltaRl=0;
	int    BdeltaR=180;
	double RdeltaPhih=3.142;
	double RdeltaPhil=0;
	int    BdeltaPhi=180;


	plot_cut_tri_mll           = new TH1D("plot_cut_tri_mll",           "plot_cut_tri_mll",          Bmll, Rmlll, Rmllh);
	plot_cut_tri_deltaR        = new TH1D("plot_cut_tri_deltaR",        "plot_cut_tri_deltaR",       BdeltaR  ,  RdeltaRl  ,  RdeltaRh);
	plot_cut_tri_deltaPhi      = new TH1D("plot_cut_tri_deltaPhi",      "plot_cut_tri_deltaPhi",     BdeltaPhi,  RdeltaPhil,  RdeltaPhih);


	plot_cut_tri_n0_nsfos_mll      = new TH1D("plot_cut_tri_n0_nsfos_mll"      ,"plot_cut_tri_n0_nsfos_mll"     , Bmll     , Rmlll     , Rmllh);
	plot_cut_tri_n0_nsfos_deltaR   = new TH1D("plot_cut_tri_n0_nsfos_deltaR"   ,"plot_cut_tri_n0_nsfos_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
	plot_cut_tri_n0_nsfos_deltaPhi = new TH1D("plot_cut_tri_n0_nsfos_deltaPhi" ,"plot_cut_tri_n0_nsfos_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

	plot_cut_tri_n0_nj_mll         = new TH1D("plot_cut_tri_n0_nj_mll"         ,"plot_cut_tri_n0_nj_mll"        , Bmll        , Rmlll        , Rmllh);
	plot_cut_tri_n0_nj_deltaR      = new TH1D("plot_cut_tri_n0_nj_deltaR"      ,"plot_cut_tri_n0_nj_deltaR"     , BdeltaR     , RdeltaRl     , RdeltaRh);
	plot_cut_tri_n0_nj_deltaPhi    = new TH1D("plot_cut_tri_n0_nj_deltaPhi"    ,"plot_cut_tri_n0_nj_deltaPhi"   , BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

	plot_cut_tri_n0_nb_mll         = new TH1D("plot_cut_tri_n0_nb_mll"         ,"plot_cut_tri_n0_nb_mll"        , Bmll        , Rmlll        , Rmllh);
	plot_cut_tri_n0_nb_deltaR      = new TH1D("plot_cut_tri_n0_nb_deltaR"      ,"plot_cut_tri_n0_nb_deltaR"     , BdeltaR     , RdeltaRl     , RdeltaRh);
	plot_cut_tri_n0_nb_deltaPhi    = new TH1D("plot_cut_tri_n0_nb_deltaPhi"    ,"plot_cut_tri_n0_nb_deltaPhi"   , BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

	plot_cut_tri_n0_dphi_mll       = new TH1D("plot_cut_tri_n0_dphi_mll"       ,"plot_cut_tri_n0_dphi_mll"      , Bmll        , Rmlll        , Rmllh);
	plot_cut_tri_n0_dphi_deltaR    = new TH1D("plot_cut_tri_n0_dphi_deltaR"    ,"plot_cut_tri_n0_dphi_deltaR"   , BdeltaR     , RdeltaRl     , RdeltaRh);
	plot_cut_tri_n0_dphi_deltaPhi  = new TH1D("plot_cut_tri_n0_dphi_deltaPhi"  ,"plot_cut_tri_n0_dphi_deltaPhi" , BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

	plot_cut_tri_n0_metpt_mll      = new TH1D("plot_cut_tri_n0_metpt_mll"      ,"plot_cut_tri_n0_metpt_mll"     , Bmll        , Rmlll        , Rmllh);
	plot_cut_tri_n0_metpt_deltaR   = new TH1D("plot_cut_tri_n0_metpt_deltaR"   ,"plot_cut_tri_n0_metpt_deltaR"  , BdeltaR     , RdeltaRl     , RdeltaRh);
	plot_cut_tri_n0_metpt_deltaPhi = new TH1D("plot_cut_tri_n0_metpt_deltaPhi" ,"plot_cut_tri_n0_metpt_deltaPhi", BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

	plot_cut_tri_n0_mll3l_mll      = new TH1D("plot_cut_tri_n0_mll3l_mll"      ,"plot_cut_tri_n0_mll3l_mll"     , Bmll        , Rmlll        , Rmllh);
	plot_cut_tri_n0_mll3l_deltaR   = new TH1D("plot_cut_tri_n0_mll3l_deltaR"   ,"plot_cut_tri_n0_mll3l_deltaR"  , BdeltaR     , RdeltaRl     , RdeltaRh);
	plot_cut_tri_n0_mll3l_deltaPhi = new TH1D("plot_cut_tri_n0_mll3l_deltaPhi" ,"plot_cut_tri_n0_mll3l_deltaPhi", BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

        plot_cut_tri_n0_m3l_mll      = new TH1D("plot_cut_tri_n0_m3l_mll"      ,"plot_cut_tri_n0_m3l_mll"     , Bmll        , Rmlll        , Rmllh);
        plot_cut_tri_n0_m3l_deltaR   = new TH1D("plot_cut_tri_n0_m3l_deltaR"   ,"plot_cut_tri_n0_m3l_deltaR"  , BdeltaR     , RdeltaRl     , RdeltaRh);
        plot_cut_tri_n0_m3l_deltaPhi = new TH1D("plot_cut_tri_n0_m3l_deltaPhi" ,"plot_cut_tri_n0_m3l_deltaPhi", BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

        plot_cut_tri_n0_mee3l_mll      = new TH1D("plot_cut_tri_n0_mee3l_mll"      ,"plot_cut_tri_n0_mee3l_mll"     , Bmll        , Rmlll        , Rmllh);
        plot_cut_tri_n0_mee3l_deltaR   = new TH1D("plot_cut_tri_n0_mee3l_deltaR"   ,"plot_cut_tri_n0_mee3l_deltaR"  , BdeltaR     , RdeltaRl     , RdeltaRh);
        plot_cut_tri_n0_mee3l_deltaPhi = new TH1D("plot_cut_tri_n0_mee3l_deltaPhi" ,"plot_cut_tri_n0_mee3l_deltaPhi", BdeltaPhi   , RdeltaPhil   , RdeltaPhih);

        plot_cut_tri_n0_mtmax_mll      = new TH1D("plot_cut_tri_n0_mtmax_mll"      ,"plot_cut_tri_n0_mtmax_mll"     , Bmll        , Rmlll        , Rmllh);
        plot_cut_tri_n0_mtmax_deltaR   = new TH1D("plot_cut_tri_n0_mtmax_deltaR"   ,"plot_cut_tri_n0_mtmax_deltaR"  , BdeltaR     , RdeltaRl     , RdeltaRh);
        plot_cut_tri_n0_mtmax_deltaPhi = new TH1D("plot_cut_tri_n0_mtmax_deltaPhi" ,"plot_cut_tri_n0_mtmax_deltaPhi", BdeltaPhi   , RdeltaPhil   , RdeltaPhih);


        plot_cut_tri_n1_nsfos_mll      = new TH1D("plot_cut_tri_n1_nsfos_mll"      ,"plot_cut_tri_n1_nsfos_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_nsfos_deltaR   = new TH1D("plot_cut_tri_n1_nsfos_deltaR"   ,"plot_cut_tri_n1_nsfos_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_nsfos_deltaPhi = new TH1D("plot_cut_tri_n1_nsfos_deltaPhi" ,"plot_cut_tri_n1_nsfos_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_nj_mll      = new TH1D("plot_cut_tri_n1_nj_mll"      ,"plot_cut_tri_n1_nj_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_nj_deltaR   = new TH1D("plot_cut_tri_n1_nj_deltaR"   ,"plot_cut_tri_n1_nj_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_nj_deltaPhi = new TH1D("plot_cut_tri_n1_nj_deltaPhi" ,"plot_cut_tri_n1_nj_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_nb_mll      = new TH1D("plot_cut_tri_n1_nb_mll"      ,"plot_cut_tri_n1_nb_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_nb_deltaR   = new TH1D("plot_cut_tri_n1_nb_deltaR"   ,"plot_cut_tri_n1_nb_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_nb_deltaPhi = new TH1D("plot_cut_tri_n1_nb_deltaPhi" ,"plot_cut_tri_n1_nb_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_pt3l_mll      = new TH1D("plot_cut_tri_n1_pt3l_mll"      ,"plot_cut_tri_n1_pt3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_pt3l_deltaR   = new TH1D("plot_cut_tri_n1_pt3l_deltaR"   ,"plot_cut_tri_n1_pt3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_pt3l_deltaPhi = new TH1D("plot_cut_tri_n1_pt3l_deltaPhi" ,"plot_cut_tri_n1_pt3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_dphi_mll      = new TH1D("plot_cut_tri_n1_dphi_mll"      ,"plot_cut_tri_n1_dphi_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_dphi_deltaR   = new TH1D("plot_cut_tri_n1_dphi_deltaR"   ,"plot_cut_tri_n1_dphi_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_dphi_deltaPhi = new TH1D("plot_cut_tri_n1_dphi_deltaPhi" ,"plot_cut_tri_n1_dphi_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_metpt_mll      = new TH1D("plot_cut_tri_n1_metpt_mll"      ,"plot_cut_tri_n1_metpt_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_metpt_deltaR   = new TH1D("plot_cut_tri_n1_metpt_deltaR"   ,"plot_cut_tri_n1_metpt_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_metpt_deltaPhi = new TH1D("plot_cut_tri_n1_metpt_deltaPhi" ,"plot_cut_tri_n1_metpt_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_mll3l_mll      = new TH1D("plot_cut_tri_n1_mll3l_mll"      ,"plot_cut_tri_n1_mll3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_mll3l_deltaR   = new TH1D("plot_cut_tri_n1_mll3l_deltaR"   ,"plot_cut_tri_n1_mll3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_mll3l_deltaPhi = new TH1D("plot_cut_tri_n1_mll3l_deltaPhi" ,"plot_cut_tri_n1_mll3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_m3l_mll      = new TH1D("plot_cut_tri_n1_m3l_mll"      ,"plot_cut_tri_n1_m3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_m3l_deltaR   = new TH1D("plot_cut_tri_n1_m3l_deltaR"   ,"plot_cut_tri_n1_m3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_m3l_deltaPhi = new TH1D("plot_cut_tri_n1_m3l_deltaPhi" ,"plot_cut_tri_n1_m3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_nsfosinz_mll      = new TH1D("plot_cut_tri_n1_nsfosinz_mll"      ,"plot_cut_tri_n1_nsfosinz_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_nsfosinz_deltaR   = new TH1D("plot_cut_tri_n1_nsfosinz_deltaR"   ,"plot_cut_tri_n1_nsfosinz_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_nsfosinz_deltaPhi = new TH1D("plot_cut_tri_n1_nsfosinz_deltaPhi" ,"plot_cut_tri_n1_nsfosinz_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n1_mt3rd_mll      = new TH1D("plot_cut_tri_n1_mt3rd_mll"      ,"plot_cut_tri_n1_mt3rd_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n1_mt3rd_deltaR   = new TH1D("plot_cut_tri_n1_mt3rd_deltaR"   ,"plot_cut_tri_n1_mt3rd_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n1_mt3rd_deltaPhi = new TH1D("plot_cut_tri_n1_mt3rd_deltaPhi" ,"plot_cut_tri_n1_mt3rd_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_nsfos_mll      = new TH1D("plot_cut_tri_n2_nsfos_mll"      ,"plot_cut_tri_n2_nsfos_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_nsfos_deltaR   = new TH1D("plot_cut_tri_n2_nsfos_deltaR"   ,"plot_cut_tri_n2_nsfos_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_nsfos_deltaPhi = new TH1D("plot_cut_tri_n2_nsfos_deltaPhi" ,"plot_cut_tri_n2_nsfos_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_nj_mll      = new TH1D("plot_cut_tri_n2_nj_mll"      ,"plot_cut_tri_n2_nj_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_nj_deltaR   = new TH1D("plot_cut_tri_n2_nj_deltaR"   ,"plot_cut_tri_n2_nj_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_nj_deltaPhi = new TH1D("plot_cut_tri_n2_nj_deltaPhi" ,"plot_cut_tri_n2_nj_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_nb_mll      = new TH1D("plot_cut_tri_n2_nb_mll"      ,"plot_cut_tri_n2_nb_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_nb_deltaR   = new TH1D("plot_cut_tri_n2_nb_deltaR"   ,"plot_cut_tri_n2_nb_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_nb_deltaPhi = new TH1D("plot_cut_tri_n2_nb_deltaPhi" ,"plot_cut_tri_n2_nb_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_pt3l_mll      = new TH1D("plot_cut_tri_n2_pt3l_mll"      ,"plot_cut_tri_n2_pt3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_pt3l_deltaR   = new TH1D("plot_cut_tri_n2_pt3l_deltaR"   ,"plot_cut_tri_n2_pt3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_pt3l_deltaPhi = new TH1D("plot_cut_tri_n2_pt3l_deltaPhi" ,"plot_cut_tri_n2_pt3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_dphi_mll      = new TH1D("plot_cut_tri_n2_dphi_mll"      ,"plot_cut_tri_n2_dphi_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_dphi_deltaR   = new TH1D("plot_cut_tri_n2_dphi_deltaR"   ,"plot_cut_tri_n2_dphi_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_dphi_deltaPhi = new TH1D("plot_cut_tri_n2_dphi_deltaPhi" ,"plot_cut_tri_n2_dphi_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_metpt_mll      = new TH1D("plot_cut_tri_n2_metpt_mll"      ,"plot_cut_tri_n2_metpt_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_metpt_deltaR   = new TH1D("plot_cut_tri_n2_metpt_deltaR"   ,"plot_cut_tri_n2_metpt_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_metpt_deltaPhi = new TH1D("plot_cut_tri_n2_metpt_deltaPhi" ,"plot_cut_tri_n2_metpt_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_mll3l_mll      = new TH1D("plot_cut_tri_n2_mll3l_mll"      ,"plot_cut_tri_n2_mll3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_mll3l_deltaR   = new TH1D("plot_cut_tri_n2_mll3l_deltaR"   ,"plot_cut_tri_n2_mll3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_mll3l_deltaPhi = new TH1D("plot_cut_tri_n2_mll3l_deltaPhi" ,"plot_cut_tri_n2_mll3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_m3l_mll      = new TH1D("plot_cut_tri_n2_m3l_mll"      ,"plot_cut_tri_n2_m3l_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_m3l_deltaR   = new TH1D("plot_cut_tri_n2_m3l_deltaR"   ,"plot_cut_tri_n2_m3l_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_m3l_deltaPhi = new TH1D("plot_cut_tri_n2_m3l_deltaPhi" ,"plot_cut_tri_n2_m3l_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);

        plot_cut_tri_n2_nsfosinz_mll      = new TH1D("plot_cut_tri_n2_nsfosinz_mll"      ,"plot_cut_tri_n2_nsfosinz_mll"     , Bmll     , Rmlll     , Rmllh);
        plot_cut_tri_n2_nsfosinz_deltaR   = new TH1D("plot_cut_tri_n2_nsfosinz_deltaR"   ,"plot_cut_tri_n2_nsfosinz_deltaR"  , BdeltaR  , RdeltaRl  , RdeltaRh);
        plot_cut_tri_n2_nsfosinz_deltaPhi = new TH1D("plot_cut_tri_n2_nsfosinz_deltaPhi" ,"plot_cut_tri_n2_nsfosinz_deltaPhi", BdeltaPhi, RdeltaPhil, RdeltaPhih);




        plot_cut_tri_mll->Sumw2();
        plot_cut_tri_deltaR->Sumw2();
        plot_cut_tri_deltaPhi->Sumw2();

        plot_cut_tri_n0_nsfos_mll->Sumw2();
        plot_cut_tri_n0_nsfos_deltaR->Sumw2();
        plot_cut_tri_n0_nsfos_deltaPhi->Sumw2();

        plot_cut_tri_n0_nj_mll->Sumw2();
        plot_cut_tri_n0_nj_deltaR->Sumw2();
        plot_cut_tri_n0_nj_deltaPhi->Sumw2();

        plot_cut_tri_n0_nb_mll->Sumw2();
        plot_cut_tri_n0_nb_deltaR->Sumw2();
        plot_cut_tri_n0_nb_deltaPhi->Sumw2();

        plot_cut_tri_n0_dphi_mll->Sumw2();
        plot_cut_tri_n0_dphi_deltaR->Sumw2();
        plot_cut_tri_n0_dphi_deltaPhi->Sumw2();

        plot_cut_tri_n0_metpt_mll->Sumw2();
        plot_cut_tri_n0_metpt_deltaR->Sumw2();
        plot_cut_tri_n0_metpt_deltaPhi->Sumw2();

        plot_cut_tri_n0_mll3l_mll->Sumw2();
        plot_cut_tri_n0_mll3l_deltaR->Sumw2();
        plot_cut_tri_n0_mll3l_deltaPhi->Sumw2();

        plot_cut_tri_n0_m3l_mll->Sumw2();
        plot_cut_tri_n0_m3l_deltaR->Sumw2();
        plot_cut_tri_n0_m3l_deltaPhi->Sumw2();

        plot_cut_tri_n0_mee3l_mll->Sumw2();
        plot_cut_tri_n0_mee3l_deltaR->Sumw2();
        plot_cut_tri_n0_mee3l_deltaPhi->Sumw2();

        plot_cut_tri_n0_mtmax_mll->Sumw2();
        plot_cut_tri_n0_mtmax_deltaR->Sumw2();
        plot_cut_tri_n0_mtmax_deltaPhi->Sumw2();


        plot_cut_tri_n1_nsfos_mll->Sumw2();
        plot_cut_tri_n1_nsfos_deltaR->Sumw2();
        plot_cut_tri_n1_nsfos_deltaPhi->Sumw2();

        plot_cut_tri_n1_nj_mll->Sumw2();
        plot_cut_tri_n1_nj_deltaR->Sumw2();
        plot_cut_tri_n1_nj_deltaPhi->Sumw2();

        plot_cut_tri_n1_nb_mll->Sumw2();
        plot_cut_tri_n1_nb_deltaR->Sumw2();
        plot_cut_tri_n1_nb_deltaPhi->Sumw2();

        plot_cut_tri_n1_pt3l_mll->Sumw2();
        plot_cut_tri_n1_pt3l_deltaR->Sumw2();
        plot_cut_tri_n1_pt3l_deltaPhi->Sumw2();

        plot_cut_tri_n1_dphi_mll->Sumw2();
        plot_cut_tri_n1_dphi_deltaR->Sumw2();
        plot_cut_tri_n1_dphi_deltaPhi->Sumw2();

        plot_cut_tri_n1_metpt_mll->Sumw2();
        plot_cut_tri_n1_metpt_deltaR->Sumw2();
        plot_cut_tri_n1_metpt_deltaPhi->Sumw2();

        plot_cut_tri_n1_mll3l_mll->Sumw2();
        plot_cut_tri_n1_mll3l_deltaR->Sumw2();
        plot_cut_tri_n1_mll3l_deltaPhi->Sumw2();

        plot_cut_tri_n1_m3l_mll->Sumw2();
        plot_cut_tri_n1_m3l_deltaR->Sumw2();
        plot_cut_tri_n1_m3l_deltaPhi->Sumw2();

        plot_cut_tri_n1_nsfosinz_mll->Sumw2();
        plot_cut_tri_n1_nsfosinz_deltaR->Sumw2();
        plot_cut_tri_n1_nsfosinz_deltaPhi->Sumw2();

        plot_cut_tri_n1_mt3rd_mll->Sumw2();
        plot_cut_tri_n1_mt3rd_deltaR->Sumw2();
        plot_cut_tri_n1_mt3rd_deltaPhi->Sumw2();


        plot_cut_tri_n2_nsfos_mll->Sumw2();
        plot_cut_tri_n2_nsfos_deltaR->Sumw2();
        plot_cut_tri_n2_nsfos_deltaPhi->Sumw2();

        plot_cut_tri_n2_nj_mll->Sumw2();
        plot_cut_tri_n2_nj_deltaR->Sumw2();
        plot_cut_tri_n2_nj_deltaPhi->Sumw2();

        plot_cut_tri_n2_nb_mll->Sumw2();
        plot_cut_tri_n2_nb_deltaR->Sumw2();
        plot_cut_tri_n2_nb_deltaPhi->Sumw2();

        plot_cut_tri_n2_pt3l_mll->Sumw2();
        plot_cut_tri_n2_pt3l_deltaR->Sumw2();
        plot_cut_tri_n2_pt3l_deltaPhi->Sumw2();

        plot_cut_tri_n2_dphi_mll->Sumw2();
        plot_cut_tri_n2_dphi_deltaR->Sumw2();
        plot_cut_tri_n2_dphi_deltaPhi->Sumw2();

        plot_cut_tri_n2_metpt_mll->Sumw2();
        plot_cut_tri_n2_metpt_deltaR->Sumw2();
        plot_cut_tri_n2_metpt_deltaPhi->Sumw2();

        plot_cut_tri_n2_mll3l_mll->Sumw2();
        plot_cut_tri_n2_mll3l_deltaR->Sumw2();
        plot_cut_tri_n2_mll3l_deltaPhi->Sumw2();

        plot_cut_tri_n2_m3l_mll->Sumw2();
        plot_cut_tri_n2_m3l_deltaR->Sumw2();
        plot_cut_tri_n2_m3l_deltaPhi->Sumw2();

        plot_cut_tri_n2_nsfosinz_mll->Sumw2();
        plot_cut_tri_n2_nsfosinz_deltaR->Sumw2();
        plot_cut_tri_n2_nsfosinz_deltaPhi->Sumw2();




//	char nameeta[100];
//	for(int cir1 = 0 ; cir1 < 12 ; cir1 ++){
//		sprintf(nameeta,"%s%d","plot_zboson_mass_vs_muon_eta",cir1);
//		plot_zboson_mass_vs_muon_eta[cir1] = new TH1D(nameeta,nameeta,140,60,130);
//		sprintf(nameeta,"%s%d","plot_zboson_mass_vs_antimuon_eta",cir1);
//		plot_zboson_mass_vs_antimuon_eta[cir1] = new TH1D(nameeta,nameeta,140,60,130);
//		plot_zboson_mass_vs_muon_eta[cir1]->Sumw2();
//		plot_zboson_mass_vs_antimuon_eta[cir1]->Sumw2();
//	}

//        plot_ee_lep_pt      = new TH1D("plot_ee_lep_pt",    "plot_ee_lep_pt",100,0,400);
//        plot_ee_lep_pt->Sumw2();
//        plot_em_lep_pt      = new TH1D("plot_em_lep_pt",    "plot_em_lep_pt",100,0,400);
//        plot_em_lep_pt->Sumw2();
//        plot_mm_lep_pt      = new TH1D("plot_mm_lep_pt",    "plot_mm_lep_pt",100,0,400);
//        plot_mm_lep_pt->Sumw2();
//        plot_ee_lep_eta     = new TH1D("plot_ee_lep_eta",   "plot_ee_lep_eta",100,-5,5);
//        plot_ee_lep_eta->Sumw2();
//        plot_em_lep_eta     = new TH1D("plot_em_lep_eta",   "plot_em_lep_eta",100,-5,5);
//        plot_em_lep_eta->Sumw2();
//        plot_mm_lep_eta     = new TH1D("plot_mm_lep_eta",   "plot_mm_lep_eta",100,-5,5);
//        plot_mm_lep_eta->Sumw2();
//        plot_ee_lep_trk_pt  = new TH1D("plot_ee_lep_trk_pt","plot_ee_lep_trk_pt",100,0,400);
//        plot_ee_lep_trk_pt->Sumw2();
//        plot_em_lep_trk_pt  = new TH1D("plot_em_lep_trk_pt","plot_em_lep_trk_pt",100,0,400);
//        plot_em_lep_trk_pt->Sumw2();
//        plot_mm_lep_trk_pt  = new TH1D("plot_mm_lep_trk_pt","plot_mm_lep_trk_pt",100,0,400);
//        plot_mm_lep_trk_pt->Sumw2();
//        plot_ee_lep_charge  = new TH1D("plot_ee_lep_charge","plot_ee_lep_charge",10,-1.0,1);
//        plot_ee_lep_charge->Sumw2();
 //       plot_em_lep_charge  = new TH1D("plot_em_lep_charge","plot_em_lep_charge",10,-1.0,1);
//        plot_em_lep_charge->Sumw2();
//        plot_mm_lep_charge  = new TH1D("plot_mm_lep_charge","plot_mm_lep_charge",10,-1.0,1);
//        plot_mm_lep_charge->Sumw2();
/*	plot_zboson_mass_forward = new TH1D("plot_zboson_mass_forward","zboson mass forward",140,60,130);
	plot_zboson_mass_forward->Sumw2();
	plot_zboson_mass_backward = new TH1D("plot_zboson_mass_backward","zboson mass backward",140,60,130);
	plot_zboson_mass_backward->Sumw2();
	char name_seed_muon[100];
	char name_seed_antimuon[100];
	char name_forward_muon[100];
	char name_forward_antimuon[100];
	for(int iplot = 0 ; iplot < 16 ; iplot++)
	{
		sprintf(name_seed_muon, "%s%d", "plot_zboson_mass_seed_muon", iplot);
		sprintf(name_seed_antimuon, "%s%d", "plot_zboson_mass_seed_antimuon", iplot);
		plot_zboson_mass_seed_muon[iplot] = new TH1D(name_seed_muon,name_seed_muon,140,60,130);
		plot_zboson_mass_seed_antimuon[iplot] = new TH1D(name_seed_antimuon,name_seed_antimuon,140,60,130);
		plot_zboson_mass_seed_muon[iplot]->Sumw2();
		plot_zboson_mass_seed_antimuon[iplot]->Sumw2();
	}
	for(int iplot = 0; iplot<16;iplot++)
	{
		sprintf(name_forward_muon,"%s%d","plot_zboson_mass_forward_muon",iplot);
		sprintf(name_forward_antimuon,"%s%d","plot_zboson_mass_forward_antimuon",iplot);
		plot_zboson_mass_forward_muon[iplot] = new TH1D(name_forward_muon, name_forward_muon,140,60,130);
		plot_zboson_mass_forward_antimuon[iplot] = new TH1D(name_forward_antimuon, name_forward_antimuon,140,60,130);
		plot_zboson_mass_forward_muon[iplot]->Sumw2();
		plot_zboson_mass_forward_antimuon[iplot]->Sumw2();
	}*/
}

void makeHists::saveHists()
{
	hf->cd();
	hf->Write();
	hf->Close();
}
 
