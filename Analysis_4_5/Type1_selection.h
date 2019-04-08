#include <iostream>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1D.h"
 
using namespace std;

      int chi2_selec( double muon_trk_chi2){
          if(muon_trk_chi2 < 4) return true;
          else  return false; 
}

      int Pt_selec( double muon_trk_pt){
          if(muon_trk_pt > 15) return true;
          else  return false;
}

      int delta_dca_selec( double muon_trk_dca, int muon_trk_nsmt){
          if(muon_trk_nsmt == 0 && fabs(muon_trk_dca) < 0.2) 
          return true;
          if(muon_trk_nsmt > 0 && fabs(muon_trk_dca) < 0.012)
          return true;
          else return false;
}
           
      int etHalo_selec( double muon_etHalo, double muon_Pt, double instlum){
          double I_scaled;
          I_scaled = ( muon_etHalo - 0.005 * instlum ) / muon_Pt ;
          if(I_scaled < 0.4 ) return true;
          else return false;
}
  
      int etTrkCone5_selec( double muon_etTrkCone5, double muon_Pt ){
          double I_scaled;
          I_scaled = muon_etTrkCone5 / muon_Pt;
          if(I_scaled  <0.4) return true;
          else return false;
}
 
     int nseg_selec( int muon_mtc_nseg ){
          if(muon_mtc_nseg > 0) return true;
          else return false;
}
    
     int trk_match_selec( int muon_trk_match ){
         if( muon_trk_match == 1 ) return true;
         else return false;
}



