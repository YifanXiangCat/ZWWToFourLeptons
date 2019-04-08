#define zeeAna_cxx
#include "tool_programs.h"
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
#include "AnglesUtil.hpp"
#include "Angle.h"
#include "global_selection.h"
#include "Type1_selection.h"


using namespace std;

void Analysis::Loop(const char* TypeName)
{

	if(isRead == false)
	{
		cout<<"this file does not content any events, skip to the next"<<endl;
		return;
	}
	
	
	bool test_on_raw = false;

	bool isData = false ;
	bool isMC = false;
	
	int debug;
	int read;

	int test =0;

	int muon1,muon2;
	double collin_cos_theta;
	double collin_phi;

	bool selec_ee;
	bool selec_em;
	bool selec_mm;
	bool selec_tri_n0;
	bool selec_tri_n1;
	bool selec_tri_n2;


	TString InputName = TString(TypeName);
	if(InputName == "data")    	isData = true;
	else if(InputName == "MC")  isMC = true;
	else	cout<<"input file type wrong! please input 'data' or 'MC'!"<<endl;
	
	if(fChain == 0)	return;
	
	Int_t Nentries = fChain->GetEntries();
	ifstream infile("read.txt");
	ofstream outfile;
	outfile.open("test.txt",ios::out);


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


	TLorentzVector lep[3];
	TLorentzVector Higgs;

	/*float lead_mu_recoid_sf;
	float subl_mu_recoid_sf;
	float lead_el_recoid_sf;
	float subl_el_recoid_sf;
	float lead_el_mva_sf;
	float subl_el_mva_sf;
	float emu_mu_recoid_sf;
	float emu_el_recoid_sf;
	float emu_el_mva_sf;
	float lead_mu_recoid_3l_sf;
	float subl_mu_recoid_3l_sf;
	float lead_el_recoid_3l_sf;
	float subl_el_recoid_3l_sf;
	float lead_el_mva_3l_sf;
	float subl_el_mva_3l_sf;
	float tert_mu_recoid_3l_sf;
	float tert_el_recoid_3l_sf;
	float tert_el_mva_3l_sf;
	float lead_mu_isoip_sf;
	float subl_mu_isoip_sf;
	float lead_el_isoip_sf;
	float subl_el_isoip_sf;
	float emu_mu_isoip_sf;
	float emu_el_isoip_sf;
	float lead_mu_isoip_3l_sf;
	float subl_mu_isoip_3l_sf;
	float lead_el_isoip_3l_sf;
	float subl_el_isoip_3l_sf;
	float tert_mu_isoip_3l_sf;
	float tert_el_isoip_3l_sf;*/

//	double trig_sf = 1;

	TRandom3 *myR = new TRandom3(13579);
	TRandom3 *MC_solenoid = new TRandom3(272727);
	cout<<"there are "<<Nentries<<" events in Loop"<<endl;
	for(int ii = 0 ; ii < Nentries ; ii++)
	{
		cout<<"pass0"<<endl;
		fChain->GetEntry(ii,0);
		cout<<"pass1"<<endl;

		int test;
		weight = 1;
		int binnumber;
		cout<<"pass2"<<endl;
		if(isMC){
			binnumber = puw_2017_central->FindBin(nTrueInt);
			weight_2017get = (double )(puw_2017_central->GetBinContent(nTrueInt));
			weight = weight * evt_scale1fb;
			weight = weight * lumi_weight;
			weight = weight * weight_2017get;
//			weight = weight * 0.002118 / 0.0047318;
//			cout<<evt_scale1fb<<"   "<<lumi_weight<<"   "<<weight_2017get<<endl;
//			cout<<"weight from puw: "<<weight_2017get<<"   scale1fb: "<<evt_scale1fb<<"   lumi: "<<lumi_weight<<"   total: "<<weight<<endl;
			
			double findpt[3];
			double findeta[3];
//			cout<<"pt num "<<(*lep_pt).size()<<"  pdg num "<<(*lep_pdgId).size()<<endl;
			for(int ii = 0 ; ii < (*lep_pt).size() ; ii ++){
				weight_lep[ii] = 0;
				if((*lep_pt)[ii]>120) findpt[ii] = 119.99;
				else findpt[ii]= (*lep_pt)[ii];
				findeta[ii] = fabs((*lep_eta)[ii]);
//				cout<<findpt[ii]<<"   "<<(*lep_pt)[ii]<<"   "<<findeta[ii]<<"   "<<(*lep_eta)[ii]<<endl;
				binnumber = SFP_muon_recoid->FindBin(findpt[ii],findeta[ii]);
				SF_muon_recoid[ii] = SFP_muon_recoid->GetBinContent(binnumber);
				binnumber = SFP_elec_recoid->FindBin(findeta[ii],findpt[ii]);
				SF_elec_recoid[ii] = SFP_elec_recoid->GetBinContent(binnumber);
				binnumber = SFP_elec_mva->FindBin(findeta[ii],findpt[ii]);//should be (*lep_eta!!!!)
				SF_elec_mva[ii] = SFP_elec_mva->GetBinContent(binnumber);
				binnumber = SFP_elec_mva_3l->FindBin(findeta[ii],findpt[ii]);//should be (*lep_eta!!!)
				SF_elec_mva_3l[ii] = SFP_elec_mva_3l->GetBinContent(binnumber);
				binnumber = SFP_muon_isoip->FindBin(findeta[ii],findpt[ii]);
				SF_muon_isoip[ii] = SFP_muon_isoip->GetBinContent(binnumber);
				binnumber = SFP_muon_isoip_3l->FindBin(findeta[ii],findpt[ii]);
				SF_muon_isoip_3l[ii] = SFP_muon_isoip_3l->GetBinContent(binnumber);
				binnumber = SFP_elec_isoip->FindBin(findeta[ii],findpt[ii]);
				SF_elec_isoip[ii] = SFP_elec_isoip->GetBinContent(binnumber);
				binnumber = SFP_elec_isoip_3l->FindBin(findeta[ii],findpt[ii]);
				SF_elec_isoip_3l[ii] = SFP_elec_isoip_3l->GetBinContent(binnumber);
//				cout<<"lep num "<<(*lep_pt).size()<<"  lep pdg"<<((*lep_pdgId)[ii])<<endl;
				if((*lep_id).at(ii) == 11) weight_lep[ii] = SF_elec_recoid[ii] * SF_elec_mva_3l[ii] * SF_elec_isoip_3l[ii];
				else weight_lep[ii] = SF_muon_recoid[ii] * SF_muon_isoip_3l[ii];
			}
			
		}

//		weight = 1;

//***************************************************//
//**                                               **//
//**   should I add presel and trigger selection?? **//
//**                                               **//
//***************************************************//
		//selection for presel*************************
//		bool presel;
//		presel = true;
//		if(firstgoodvertex != 0)       presel = false;
//		if(Flag_AllEventFilters <= 0)  presel = false;
//		if(vetophoton != 0)            presel = false;
//		if(evt_passgoodrunlist <= 0 )  presel = false;
//		if(nVlep < 2)                  presel = false;
//		if(nLlep < 2)                  presel = false;

//		if(presel == true)   count_presel +=weight;

//		if(presel == false)  continue;


		//selection for trigger!!!******************************
//		bool passtrigger;
//		passtrigger = true;
//		if(passTrigger == 0)              passtrigger = false;
//		if(pass_duplicate_ee_em_mm == 0)  passtrigger = false;
		
//		if(passtrigger == true)   count_trig +=weight;
//		if(passtrigger == false) continue;

		//selection for isdilepton******************

		//selection for istrilepton***************

		//selection for isfourlepton***************
		int number_lep = (*lep_pt).size();
		cout<<number_lep<<endl;



		//cout<<"test = "<<test<<"   "<<MTmax<<endl;
	}//end of loop

}//end of whole function
































