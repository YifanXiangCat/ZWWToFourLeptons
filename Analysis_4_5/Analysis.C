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

	int count_0lep=0;
	int count_1lep=0;
	int count_2lep=0;
	int count_3lep=0;
	int count_4lep=0;
	int count_5lep=0;


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
		fChain->GetEntry(ii,0);

		int test;
		weight = 1;
		int binnumber;
/*		if(isMC){
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
*/
//		weight = 1;

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
		int  number_lep = (*lep_pt).size();
		if(number_lep == 0 ) count_0lep ++;
		else if(number_lep == 1) count_1lep++;
		else if(number_lep == 2) count_2lep++;
		else if(number_lep == 3) count_3lep++;
		else if(number_lep == 4) count_4lep++;
		else if(number_lep == 5) count_5lep++;


		if(number_lep < 4 )  continue;

		bool passfourleptons;
		passfourleptons = true;
		TLorentzVector leptons[10];
		TLorentzVector dilep;
		for(int jj = 0 ; jj < number_lep ; jj ++){
			leptons[jj].SetPtEtaPhiE(  ( (*lep_pt).at(jj) ),( (*lep_eta).at(jj) ),( (*lep_phi).at(jj) ),( (*lep_energy).at(jj) )  );
		}
		int lepz1,lepz2;
		lepz1 = -999;
		lepz2 = -990;
		double compare;
		compare = 10 ;
		bool ifpass;
		int id1;
		int id2;
		bool ifpassselection;
		ifpassselection = true;
		//find the first two particles decay from a Z boson
		for(int jj = 0 ; jj < (number_lep-1) ; jj ++){
			for(int kk = jj + 1 ; kk < number_lep ; kk++){
				dilep = leptons[jj] + leptons[kk];
				id1 = (*lep_id).at(jj);
				id2 = (*lep_id).at(kk);
				ifpass = true;
				if(fabs(dilep.M()-91.1875) > 10)  ifpass = false;
				if(abs(id1)!= abs(id2)) ifpass = false;
				if(id1 == id2)  ifpass = false;
				if(ifpass && fabs(dilep.M() - 91.1875)> compare) ifpass = false;
				if( ifpass ){
					if(compare!=10)  cout<<"more than one di-lepton pair"<<endl;
					compare = fabs(dilep.M()-91.1876);
					lepz1 = jj;
					lepz2 = kk;
				}
			}
		}

		if(lepz1 == -999){//if I do not find two particle decay from a zboson
			ifpassselection = false;
		}
		bool ifsamecharge;
		//the total charge of the four leptons should be 0
		int total_charge;
		total_charge = 0;
		for(int jj = 0 ; jj<(number_lep) ; jj++){
			int id = (*lep_id).at(jj);
			int charge = id/fabs(id);
			total_charge = total_charge + charge;
		}
		if(total_charge != 0) ifpassselection = false;
		//any different charge leptons should have an invarient mass higher than 12
		for(int jj = 0 ; jj < (number_lep-1) ; jj++){
			for(int kk = jj+1 ; kk <number_lep ; kk++){
				ifpass = true;
				id1 = (*lep_id).at(jj);
				id2 = (*lep_id).at(kk);
				dilep = leptons[jj] + leptons[kk];
				if(id1 > 0 && id2 > 0) ifsamecharge = true;
				else if(id1 < 0 && id2 < 0) ifsamecharge = true;
				else ifsamecharge = false;
				if(!ifsamecharge)  
					if(dilep.M() < 12)  ifpassselection = false;
			}
		}
		//any event should have nb=0
		if(nb!=0) ifpassselection = false;
		if(ifpassselection == false){
		       continue;
		}



		//cout<<"test = "<<test<<"   "<<MTmax<<endl;
	}//end of loop

	cout<<"0 lep number : "<<count_0lep<<endl;
	cout<<"1 lep number : "<<count_1lep<<endl;
	cout<<"2 lep number : "<<count_2lep<<endl;
	cout<<"3 lep number : "<<count_3lep<<endl;
	cout<<"4 lep number : "<<count_4lep<<endl;
	cout<<"5 lep number : "<<count_5lep<<endl;

}//end of whole function
































