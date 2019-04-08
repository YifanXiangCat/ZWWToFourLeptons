#include <iostream>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1D.h"

using namespace std;

      int charge_selec( int muon_trk_charge1 , int muon_trk_charge2 ){
          if( muon_trk_charge1 == 1 && muon_trk_charge2 == -1) return 1;
          else if( muon_trk_charge1 == -1 && muon_trk_charge2 == 1) return 2;
          else return false;
}
      
      int trk_dz_selec( double muon_trk_z1, double muon_trk_z2 , 
                        int muon_trk_nsmt1, int muon_trk_nsmt2){
          double trk_dz;
          int N_smt;
      
          if(muon_trk_nsmt1 >=2 && muon_trk_nsmt2 >=2) N_smt = 2;      // get the N_smt
          else if(muon_trk_nsmt1 < 2 && muon_trk_nsmt2 < 2) N_smt = 0;   
          else N_smt = 1;
   
          trk_dz = fabs( muon_trk_z1 - muon_trk_z2);
          
          if( N_smt == 0 && trk_dz < 2) return true;
          if( N_smt == 1 && trk_dz < 1.5) return true;
          if( N_smt == 2 && trk_dz < 1) return true;
          return false;
}
      int cos_theta_selec( double muon_trk_px1, double muon_trk_py1, double muon_trk_pz1,
                           double muon_trk_px2, double muon_trk_py2, double muon_trk_pz2){
          double cos_theta = 0;
          double p1_p2, p1, p2 = 0;
          p1_p2 = muon_trk_px1 * muon_trk_px2 + muon_trk_py1 * muon_trk_py2 + muon_trk_pz1 * muon_trk_pz2 ;
          p1 = sqrt( muon_trk_px1 * muon_trk_px1 + muon_trk_py1 * muon_trk_py1 + muon_trk_pz1 * muon_trk_pz1 );
          p2 = sqrt( muon_trk_px2 * muon_trk_px2 + muon_trk_py2 * muon_trk_py2 + muon_trk_pz2 * muon_trk_pz2 );
          cos_theta = p1_p2 / (p1 * p2);
          
          if( cos_theta > -0.99985) return true;
          else return false;
}






