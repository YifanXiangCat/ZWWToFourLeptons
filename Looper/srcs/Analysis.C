#define zeeAna_cxx
#include "Analysis.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
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
#include <TGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <iomanip>


using namespace std;

void Analysis::Loop(const char* TypeName)
{

    if (isRead == false)
    {
        cout << "this file does not content any events, skip to the next" << endl;
        return;
    }


    bool test_on_raw = false;

    bool isData = false ;
    bool isMC = false;

    int debug;
    int read;

    int test = 0;

    int muon1, muon2;
    double collin_cos_theta;
    double collin_phi;

    bool selec_ee;
    bool selec_em;
    bool selec_mm;
    bool selec_tri_n0;
    bool selec_tri_n1;
    bool selec_tri_n2;

    double count_0lep = 0;
    double count_1lep = 0;
    double count_2lep = 0;
    double count_3lep = 0;
    double count_4lep = 0;
    double count_5lep = 0;


    TString InputName = TString(TypeName);
    if (InputName == "data")    	isData = true;
    else if (InputName == "MC")  isMC = true;
    else	cout << "input file type wrong! please input 'data' or 'MC'!" << endl;

    if (fChain == 0)	return;

    Int_t Nentries = fChain->GetEntries();
    ifstream infile("read.txt");
    ofstream outfile;
    outfile.open("test.txt", ios::out);


    double weight = 0;
    double weight_lep[3] ;
    double weight_2017get = 0;

    float SF_muon_recoid[3];
    float SF_elec_recoid[3];
    float SF_elec_mva[3];
    float SF_elec_mva_3l[3];
    float SF_muon_isoip[3];
    float SF_muon_isoip_3l[3];
    float SF_elec_isoip[3];
    float SF_elec_isoip_3l[3];

    double count_notpassselection = 0;
    double count_passselection = 0;
    double count_morethanone = 0;
    double count_one_Z_lep_fake = 0;
    double count_more_than_5_lep = 0;
    TRandom3 *myR = new TRandom3(13579);
    TRandom3 *MC_solenoid = new TRandom3(272727);
    cout << "there are " << Nentries << " events in Loop" << endl;
    for (int ii = 0 ; ii < Nentries ; ii++)
    {
        fChain->GetEntry(ii, 0);
        int test;
        weight = 1;
        int binnumber;
        weight = weight * evt_scale1fb * 137;
        //selection for isfourlepton***************
        int  number_lep = (*lep_pt).size();
        if (number_lep == 0) count_0lep += weight;
        else if (number_lep == 1) count_1lep += weight;
        else if (number_lep == 2) count_2lep += weight;
        else if (number_lep == 3) count_3lep += weight;
        else if (number_lep == 4) count_4lep += weight;
        else if (number_lep == 5) count_5lep + weight;
        TLorentzVector leptons[10];
        TLorentzVector dilep;
        TLorentzVector dilep34;
        for (int jj = 0 ; jj < number_lep ; jj ++)
            leptons[jj].SetPtEtaPhiE(((*lep_pt).at(jj)), ((*lep_eta).at(jj)), ((*lep_phi).at(jj)), ((*lep_energy).at(jj)));
        int lepz1, lepz2, lep3, lep4, lep5;
        lepz1 = -999;
        lepz2 = -999;
        lep3 = -999;
        lep4 = -999;
        lep5 = -999;
        double compare;
        compare = 10 ;
        bool ifpass;
        int id1;
        int id2;
        bool ifpassselection;
        ifpassselection = true;
        //find the first two particles decay from a Z boson
        for (int jj = 0 ; jj < (number_lep - 1) ; jj ++)
        {
            for (int kk = jj + 1 ; kk < number_lep ; kk++)
            {
                dilep = leptons[jj] + leptons[kk];
                id1 = (*lep_id).at(jj);
                id2 = (*lep_id).at(kk);
                ifpass = true;
                if (fabs(dilep.M() - 91.1875) > 10)  ifpass = false;
                if (abs(id1) != abs(id2)) ifpass = false;
                if (id1 == id2)  ifpass = false;
                if (ifpass && fabs(dilep.M() - 91.1875) > compare) ifpass = false;
                if (ifpass)
                {
                    if (compare != 10)  count_morethanone += weight;
                    compare = fabs(dilep.M() - 91.1876);
                    lepz1 = jj;
                    lepz2 = kk;
                }
            }
        }
        if (lepz1 == -999)//if I do not find two particle decay from a zboson
            ifpassselection = false;
        if (ifpassselection == false)  continue;
        dilep = leptons[lepz1] + leptons[lepz2];
        bool iffake[2];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
                myhists->plot_dilep_good_Mll[0]->Fill(dilep.M(), weight);
            else if (iffake[0] || iffake[1])
            {
                myhists->plot_dilep_fake_Mll[0]->Fill(dilep.M(), weight);
                if ((*lep_mc_id).at(lepz1) == 1 && (*lep_mc_id).at(lepz2) == 2)   cout << (*lep_mc_id).at(lepz1) << "   " << (*lep_mc_id).at(lepz2) << endl;
                if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 1)   cout << (*lep_mc_id).at(lepz1) << "   " << (*lep_mc_id).at(lepz2) << endl;
            }
            else
            {
                myhists->plot_dilep_rest_Mll[0]->Fill(dilep.M(), weight);
            }
        }
        bool ifsamecharge;
        int count_pass_relIso ;
        count_pass_relIso = 0;
        //ask all the leptosn relIso < 0.1
        for (int jj = 0 ; jj < number_lep ; jj++)
        {
            if ((*lep_relIso03EA).at(jj) > 0.4)
            {
                if (jj == lepz1)   count_one_Z_lep_fake += weight;
                if (jj == lepz2)   count_one_Z_lep_fake += weight;
                continue;
            }
            if (jj == lepz1)
            {
                count_pass_relIso ++ ;
                continue;
            }
            if (jj == lepz2)
            {
                count_pass_relIso ++;
                continue;
            }
            if (lep3 == -999  && (*lep_relIso03EA).at(jj) < 0.1)
            {
                count_pass_relIso ++;
                lep3 = jj;
                continue;
            }
            if (lep4 == -999 && (*lep_relIso03EA).at(jj) < 0.1)
            {
                count_pass_relIso ++;
                lep4 = jj;
                continue;
            }
            if (lep5 == -999)
            {
                count_pass_relIso ++;
                lep5 = jj;
                count_more_than_5_lep += weight;
                ifpassselection = false;
            }
        }
        if (count_pass_relIso != 4)  ifpassselection = false; //if I don't find the 4th lepton
        dilep = leptons[lepz1] + leptons[lepz2];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
                myhists->plot_dilep_good_Mll[1]->Fill(dilep.M(), weight);
            else if (iffake[0] || iffake[1])
                myhists->plot_dilep_fake_Mll[1]->Fill(dilep.M(), weight);
            else
                myhists->plot_dilep_rest_Mll[1]->Fill(dilep.M(), weight);
        }
        int lep_find[4];
        lep_find[0] = lepz1;
        lep_find[1] = lepz2;
        lep_find[2] = lep3;
        lep_find[3] = lep4;
        if (!ifpassselection) continue;
        //the total charge of the four leptons should be 0
        int total_charge;
        total_charge = 0;
        for (int jj = 0 ; jj < 4 ; jj ++)
        {
            int id = (*lep_id).at(lep_find[jj]);
            int charge = id / fabs(id);
            total_charge = total_charge + charge;
        }
        if (total_charge != 0) ifpassselection = false;
        dilep = leptons[lepz1] + leptons[lepz2];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
                myhists->plot_dilep_good_Mll[2]->Fill(dilep.M(), weight);
            else if (iffake[0] || iffake[1])
                myhists->plot_dilep_fake_Mll[2]->Fill(dilep.M(), weight);
            else
                myhists->plot_dilep_rest_Mll[2]->Fill(dilep.M(), weight);
        }
        if (!ifpassselection) continue;
        //the highest pt > 25 second > 20
        double HighestPt1 = 0 ;
        double HighestPt2 = 0 ;
        for (int jj = 0 ; jj < 4 ; jj++)
        {
            if (HighestPt1 < (*lep_pt).at(lep_find[jj]))
            {
                HighestPt2 = HighestPt1;
                HighestPt1 = (*lep_pt).at(lep_find[jj]);
                continue;
            }
            if (HighestPt2 < (*lep_pt).at(lep_find[jj]))      HighestPt2 = (*lep_pt).at(lep_find[jj]);
        }
        if (leptons[lepz1].Pt() > leptons[lepz2].Pt())
        {
            if (leptons[lepz1].Pt() < 25)  ifpassselection = false;
        }
        else  if (leptons[lepz2].Pt() < 25)  ifpassselection = false;
        dilep = leptons[lepz1] + leptons[lepz2];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
                myhists->plot_dilep_good_Mll[3]->Fill(dilep.M(), weight);
            else if (iffake[0] || iffake[1])
                myhists->plot_dilep_fake_Mll[3]->Fill(dilep.M(), weight);
            else
                myhists->plot_dilep_rest_Mll[3]->Fill(dilep.M(), weight);
        }
        //any different charge leptons should have an invarient mass higher than 12
        for (int jj = 0 ; jj < 3 ; jj ++)
        {
            for (int kk = jj + 1 ; kk < 4 ; kk ++)
            {
                ifpass = true;
                id1 = (*lep_id).at(lep_find[jj]);
                id2 = (*lep_id).at(lep_find[kk]);
                dilep = leptons[jj] + leptons[kk];
                if (id1 > 0 && id2 < 0) ifsamecharge = true;
                else ifsamecharge = false;
                if (!ifsamecharge)
                    if (dilep.M() < 12)   ifpassselection = false;
            }
        }
        dilep = leptons[lepz1] + leptons[lepz2];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
                myhists->plot_dilep_good_Mll[4]->Fill(dilep.M(), weight);
            else if (iffake[0] || iffake[1])
                myhists->plot_dilep_fake_Mll[4]->Fill(dilep.M(), weight);
            else
                myhists->plot_dilep_rest_Mll[4]->Fill(dilep.M(), weight);
        }
        //any event should have nb=0
        if (nb != 0) ifpassselection = false;
        dilep = leptons[lepz1] + leptons[lepz2];
        dilep34 = leptons[lep3] + leptons[lep4];
        iffake[0] = !((*lep_mc_id).at(lepz1) == 2 || (*lep_mc_id).at(lepz1) == 1);
        iffake[1] = !((*lep_mc_id).at(lepz2) == 2 || (*lep_mc_id).at(lepz2) == 1);
        if (ifpassselection)
        {
            if ((*lep_mc_id).at(lepz1) == 2 && (*lep_mc_id).at(lepz2) == 2)
            {
                myhists->plot_dilep_good_Mll[5]->Fill(dilep.M(), weight);
                myhists->plot_dilep_good_relIso->Fill((*lep_relIso03EA).at(lepz1), weight);
                myhists->plot_dilep_good_relIso->Fill((*lep_relIso03EA).at(lepz2), weight);
                if (fabs((*lep_id).at(lep3)) != fabs((*lep_id).at(lep4)))    myhists->plot_dilep_good_Mll_em->Fill(dilep34.M(), weight);
                if (fabs((*lep_id).at(lep3)) == fabs((*lep_id).at(lep4)))
                {
                    if (fabs(dilep34.M() - 91.1876) < 10)   myhists->plot_dilep_good_Mll_onZ->Fill(dilep34.M(), weight);
                    else   myhists->plot_dilep_good_Mll_offZ->Fill(dilep34.M(), weight);
                }
            }
            else if ((*lep_mc_id).at(lepz1) < 0 || (*lep_mc_id).at(lepz2) < 0)
            {
                myhists->plot_dilep_fake_Mll[5]->Fill(dilep.M(), weight);
                if ((*lep_mc_id).at(lepz1) != 2)
                    myhists->plot_dilep_fake_relIso->Fill((*lep_relIso03EA).at(lepz1), weight);
                if ((*lep_mc_id).at(lepz2) != 2)
                    myhists->plot_dilep_fake_relIso->Fill((*lep_relIso03EA).at(lepz2), weight);
                if (fabs((*lep_id).at(lep3)) != fabs((*lep_id).at(lep4)))    myhists->plot_dilep_fake_Mll_em->Fill(dilep34.M(), weight);
                if (fabs((*lep_id).at(lep3)) == fabs((*lep_id).at(lep4)))
                {
                    if (fabs(dilep34.M() - 91.1876) < 10)   myhists->plot_dilep_fake_Mll_onZ->Fill(dilep34.M(), weight);
                    else   myhists->plot_dilep_fake_Mll_offZ->Fill(dilep34.M(), weight);
                }
            }
            else
            {
                myhists->plot_dilep_rest_Mll[5]->Fill(dilep.M(), weight);
                if (fabs((*lep_id).at(lep3)) != fabs((*lep_id).at(lep4)))    myhists->plot_dilep_rest_Mll_em->Fill(dilep34.M(), weight);
                if (fabs((*lep_id).at(lep3)) == fabs((*lep_id).at(lep4)))
                {
                    if (fabs(dilep34.M() - 91.1876) < 10)   myhists->plot_dilep_rest_Mll_onZ->Fill(dilep34.M(), weight);
                    else   myhists->plot_dilep_rest_Mll_offZ->Fill(dilep34.M(), weight);
                }
            }
        }
        double St;//calculation for St
        St = ht + leptons[lepz1].Pt() + leptons[lepz2].Pt() + leptons[lep3].Pt() + leptons[lep4].Pt() + met_pt;
        TLorentzVector myjets[30];
        int  number_jet = (*jets_pt).size();
        for (int jj = 0 ; jj < number_jet ; jj ++)
            myjets[jj].SetPtEtaPhiM(((*jets_pt).at(jj)), ((*jets_eta).at(jj)), ((*jets_phi).at(jj)), ((*jets_mass).at(jj)));
        double Mjet_DeltaR ;
        int jet_DR1, jet_DR2, jet_Pt1, jet_Pt2;
        jet_DR1 = -999;
        jet_DR2 = -999;
        jet_Pt1 = 0;
        jet_Pt2 = 0;
        //find leading Pt jets;
        if (nj  > 1)
        {
            for (int jj = 0 ; jj < number_jet ; jj ++)
            {
                if (myjets[jj].Pt() < 30)  continue;
                if (myjets[jj].Eta() > 2.5)  continue;
                if (myjets[jj].Eta() < -2.5)  continue;
                if (myjets[jet_Pt1].Pt() < myjets[jj].Pt())
                {
                    jet_Pt2 = jet_Pt1;
                    jet_Pt1 = jj;
                    continue;
                }
                if (myjets[jet_Pt2].Pt() < myjets[jj].Pt())
                    jet_Pt2 = jj;
            }
        }
        //find leading DR pair jets and their mjets;
        double jet_min_DR;
        jet_min_DR = 100000;
        double jet_min_DR_Mjets = 0;
        TLorentzVector dijets;
        for (int jj = 0 ; jj < number_jet - 1 ; jj ++)
        {
            for (int kk = jj ; kk < number_jet ; kk ++)
            {
                dijets = myjets[jj] + myjets[kk];
                if (myjets[jj].DeltaR(myjets[kk]) < jet_min_DR)   jet_min_DR_Mjets = dijets.M();
            }
        }
        double jet_max_Pt_Mjets;
        dijets = myjets[jet_Pt1] + myjets[jet_Pt2];
        jet_max_Pt_Mjets = dijets.M();
        int npflep[8];
        npflep[0] = -999;
        npflep[1] = -999;
        npflep[2] = -999;
        npflep[3] = -999;
        for (int jj = 0 ; jj < 4 ; jj ++)
        {
            double comparePt;
            comparePt = leptons[jj].Pt();
            for (int kk = 0 ; kk < 4 ; kk ++)
            {
                if (npflep[kk] == -999)
                {
                    npflep[kk] = jj;
                    break;
                }
                if (leptons[npflep[kk]].Pt() < comparePt)
                {
                    for (int mm = 4 ; mm > kk ; mm --)
                        npflep[mm] = npflep[mm - 1];
                    npflep[kk] = jj;
                    break;
                }
            }
        }
        if (leptons[lepz1].Pt() > leptons[lepz2].Pt())
        {
            npflep[4] = lepz1;
            npflep[5] = lepz2;
        }
        else
        {
            npflep[4] = lepz2;
            npflep[5] = lepz1;
        }
        if (leptons[lep3].Pt() > leptons[lep4].Pt())
        {
            npflep[6] = lep3;
            npflep[7] = lep4;
        }
        else
        {
            npflep[6] = lep4;
            npflep[7] = lep3;
        }
        if (fabs((*lep_sip3d).at(lep3)) > 5) ifpassselection = false;
        if (fabs((*lep_sip3d).at(lep4)) > 5) ifpassselection = false;
        if ((*lep_relIso03EA).at(lepz1) > 0.2)  ifpassselection = false;
        if ((*lep_relIso03EA).at(lepz2) > 0.2)  ifpassselection = false;
        if (ifpassselection)
        {
            dilep34 = leptons[lep3] + leptons[lep4];
            if (fabs((*lep_id).at(lep3)) != fabs((*lep_id).at(lep4)))      //em plots
            {
                myhists->plot_lepton_Mll[0]->Fill(dilep34.M(), weight);
                myhists->plot_lepton_missingET[0]->Fill(met_pt, weight);
                myhists->plot_nj[0]->Fill(nj, weight);
                myhists->plot_ht[0]->Fill(ht, weight);
                myhists->plot_St[0]->Fill(St, weight);
                if (nj > 1)
                {
                    myhists->plot_Mjet_DeltaR[0]->Fill(jet_min_DR_Mjets, weight);
                    myhists->plot_Mjet_LeadPt[0]->Fill(jet_max_Pt_Mjets, weight);
                }
                for (int kk = 0 ; kk < 8 ; kk ++)
                {
                    myhists->plot_lepton_Pt[0][kk]->Fill((*lep_pt).at(npflep[kk]), weight);
                    myhists->plot_lepton_Eta[0][kk]->Fill((*lep_eta).at(npflep[kk]), weight);
                    myhists->plot_lepton_relIsoEA[0][kk]->Fill((*lep_relIso03EA).at(npflep[kk]), weight);
                    myhists->plot_lepton_sip3d[0][kk]->Fill((*lep_sip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_ip[0][kk]->Fill((*lep_ip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_dxy[0][kk]->Fill((*lep_dxy).at(npflep[kk]), weight);
                    myhists->plot_lepton_dz[0][kk]->Fill((*lep_dz).at(npflep[kk]), weight);
                    myhists->plot_lepton_mva[0][kk]->Fill((*lep_mva).at(npflep[kk]), weight);
                }
            }
            else if (fabs(dilep34.M() - 91.1876) < 10) // on Z plot
            {
                myhists->plot_lepton_Mll[1]->Fill(dilep34.M(), weight);
                myhists->plot_lepton_missingET[1]->Fill(met_pt, weight);
                myhists->plot_nj[1]->Fill(nj, weight);
                myhists->plot_ht[1]->Fill(ht, weight);
                myhists->plot_St[1]->Fill(St, weight);
                if (nj > 1)
                {
                    myhists->plot_Mjet_DeltaR[1]->Fill(jet_min_DR_Mjets, weight);
                    myhists->plot_Mjet_LeadPt[1]->Fill(jet_max_Pt_Mjets, weight);
                }
                for (int kk = 0 ; kk < 8 ; kk ++)
                {
                    myhists->plot_lepton_Pt[1][kk]->Fill((*lep_pt).at(npflep[kk]), weight);
                    myhists->plot_lepton_Eta[1][kk]->Fill((*lep_eta).at(npflep[kk]), weight);
                    myhists->plot_lepton_relIsoEA[1][kk]->Fill((*lep_relIso03EA).at(npflep[kk]), weight);
                    myhists->plot_lepton_sip3d[1][kk]->Fill((*lep_sip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_ip[1][kk]->Fill((*lep_ip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_dxy[1][kk]->Fill((*lep_dxy).at(npflep[kk]), weight);
                    myhists->plot_lepton_dz[1][kk]->Fill((*lep_dz).at(npflep[kk]), weight);
                    myhists->plot_lepton_mva[1][kk]->Fill((*lep_mva).at(npflep[kk]), weight);
                }
            }
            else  //off Z plots
            {
                myhists->plot_lepton_Mll[2]->Fill(dilep34.M(), weight);
                myhists->plot_lepton_missingET[2]->Fill(met_pt, weight);
                myhists->plot_nj[2]->Fill(nj, weight);
                myhists->plot_ht[2]->Fill(ht, weight);
                myhists->plot_St[2]->Fill(St, weight);
                if (nj > 1)
                {
                    myhists->plot_Mjet_DeltaR[2]->Fill(jet_min_DR_Mjets, weight);
                    myhists->plot_Mjet_LeadPt[2]->Fill(jet_max_Pt_Mjets, weight);
                }
                for (int kk = 0 ; kk < 8 ; kk ++)
                {
                    myhists->plot_lepton_Pt[2][kk]->Fill((*lep_pt).at(npflep[kk]), weight);
                    myhists->plot_lepton_Eta[2][kk]->Fill((*lep_eta).at(npflep[kk]), weight);
                    myhists->plot_lepton_relIsoEA[2][kk]->Fill((*lep_relIso03EA).at(npflep[kk]), weight);
                    myhists->plot_lepton_sip3d[2][kk]->Fill((*lep_sip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_ip[2][kk]->Fill((*lep_ip3d).at(npflep[kk]), weight);
                    myhists->plot_lepton_dxy[2][kk]->Fill((*lep_dxy).at(npflep[kk]), weight);
                    myhists->plot_lepton_dz[2][kk]->Fill((*lep_dz).at(npflep[kk]), weight);
                    myhists->plot_lepton_mva[2][kk]->Fill((*lep_mva).at(npflep[kk]), weight);
                }
            }
        }
        if (ifpassselection == false)
        {
            count_notpassselection += weight;
            continue;
        }
        count_passselection += weight;

    }//end of loop

    cout << "0 lep number : " << count_0lep << endl;
    cout << "1 lep number : " << count_1lep << endl;
    cout << "2 lep number : " << count_2lep << endl;
    cout << "3 lep number : " << count_3lep << endl;
    cout << "4 lep number : " << count_4lep << endl;
    cout << "5 lep number : " << count_5lep << endl << endl;
    cout << "pass selection : " << count_passselection << endl;
    cout << "not pass selec : " << count_notpassselection << endl;




}//end of whole function

// eof
