#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
using namespace std;

//CB function
Double_t flipCrystalBall(Double_t mean, Double_t sigma, Double_t alpha, Double_t n, TRandom* r)
{
  // function to return random numbers according to Crystal Ball function
  // lifted directly from wz_epmcs/src/WZ_Utils.hpp
  // modified by mikav:
  // if alpha or n are less than 0, just return a gaussian
  // random number
  if(alpha <= 0 || n <= 0){
    return mean + r->Gaus(0,sigma);
  }

  //from here, the function is unchanged from WZ_Utils
//  assert(alpha>0);

  // Auxiliary calculations
  double A=pow(n/fabs(alpha),n)*exp(-alpha*alpha/2.);
  double B=n/fabs(alpha)-fabs(alpha);

  // Surfaces of the two halves
  //
  // gaussian part one
  double Sgauss=sigma*sqrt(TMath::Pi()/2.);
  // gaussian part two
  Sgauss+=sigma*sqrt(TMath::Pi()/2.)*TMath::Erf(alpha/sqrt(2.));
  // tail
  double Stail=A*sigma/(n-1.)*pow(B+alpha,-n+1.);

  // Pick one of the halves
  double S=Stail+Sgauss;
  double ran=r->Rndm();
  if (ran<(Stail/S)) {
    // generate a tail event
    double d=r->Rndm();
    return sigma*B+mean-sigma*(B+alpha)*pow(d,1./(-n+1.));
  } else {
    // generate a gaussian event
    while ((ran=r->Gaus(mean,sigma))<(mean-alpha*sigma)) {
      // flip again
    }
    return ran;
  }
}

//selection programs

//eta selection
bool InCC(double eta, double eta_cuts_low, double eta_cuts_high)
{
 bool pass=false;
 if(fabs(eta)<=eta_cuts_high && fabs(eta)>= eta_cuts_low) pass=true;
 return pass;
}

bool InEC(double eta, double eta_cuts_low, double eta_cuts_high)
{
 bool pass=false;
 if(fabs(eta)<=eta_cuts_high && fabs(eta)>= eta_cuts_low) pass=true;
 return pass;
}

//pt selection
bool passPt(double pt, double pt_cuts)
{
 bool pass=false;
 if(pt>=pt_cuts) pass=true;
 return pass;
}

//fiducial selection
bool passEtaFiducial(int fiducial, int cut_fiducial)
{
 bool pass=false;
 if(cut_fiducial==1 && fiducial==1) pass=true;
 if(cut_fiducial==0) pass=true;
 return pass;
}

bool passPhiFiducial(double phimod, int type)
{
 bool pass=false;
 if(type==0) pass=true;
 else if(type==1)//select phimod 0.1~0.9
  {
   if(phimod>=0.1 && phimod<=0.9) pass=true;
  }
 else if(type==2)//select phimod 0~0.1 and 0.9~1.0
  {
   if(phimod<0.1 || phimod>0.9) pass=true;
  }
 else pass=true;
 return pass;
}

//pre-emid cuts
bool passIsoEmf(double iso, double emf, double iso_cuts, double emf_cuts)
{
 bool pass=false;
 if(iso<=iso_cuts && emf>=emf_cuts) pass=true;
 return pass;
}

//hmx cuts
bool passHmx(double hmx, double hmx_cuts)
{
 bool pass=false;
 if(hmx<=hmx_cuts) pass=true;
 if(hmx_cuts==999) pass=true;
 return pass;
}

//onn cuts
bool passOnn(double onn, double onn_cuts)
{
 bool pass=false;
 if(onn>=onn_cuts) pass=true;
 return pass;
}

//full emid cuts
bool passIsoEmfHmxOnnTrkiso(double iso, double emf, double hmx, double onn, double trkiso, double iso_cuts, double emf_cuts, double hmx_cuts, double onn_cuts, double trkiso_cuts)
{
 bool pass=false;
 if(iso<=iso_cuts && emf>=emf_cuts && hmx<=hmx_cuts && onn>=onn_cuts && trkiso<=trkiso_cuts) pass=true;
 return pass;
}

bool passIsoEmfHmxOnnTrkisoSigphi(double iso, double emf, double hmx, double onn, double trkiso, double sigphi, double iso_cuts, double emf_cuts, double hmx_cuts, double onn_cuts, double trkiso_cuts, double sigphi_cuts)
{
 bool pass=false;
 if(iso<=iso_cuts && emf>=emf_cuts && hmx<=hmx_cuts && onn>=onn_cuts && trkiso<=trkiso_cuts && sigphi<=sigphi_cuts) pass=true;
 return pass;
 
}

//trk matching
bool passTrkQuality(int type, double trkprob, double trkchi2, double trkpt, double rdca, double trksign, int nsmt, int ncft, double trkprob_cuts[], double trkchi2_cuts[], double trkpt_cuts[], double rdca_cuts[], double trksign_cuts[], int nsmt_cuts[], int ncft_cuts[], int trk_select[], double eta)
{
 bool pass=false;
 if(type==1 || type==2 || type==3 || type==4)
  {
   if(trkprob>=trkprob_cuts[type-1] && trkchi2<=trkchi2_cuts[type-1] && 
      trkpt>=trkpt_cuts[type-1] && fabs(rdca)<=rdca_cuts[type-1] && 
      fabs(trksign)>=trksign_cuts[type-1] && nsmt>=nsmt_cuts[type-1] && 
      ncft>=ncft_cuts[type-1])
    pass=true;
  }
 if(trk_select[0]==0 && fabs(eta)<=1.1) pass=true;
 if(trk_select[1]==0 && fabs(eta)>=1.5) pass=true;
 return pass;
}

//jet selection
bool passJetHmxTrkIso(double hmx, double trkiso, double hmx_cuts, double trkiso_cuts)
{
 bool pass=false;
 if(hmx>hmx_cuts && trkiso>trkiso_cuts) pass=true;
 return pass;
}

//QCD selectron
bool passQCDHmxTrkIso(double hmx, double trkiso, double hmx_cuts, double trkiso_cuts)
{
 bool pass=false;
 if(hmx>hmx_cuts || trkiso>trkiso_cuts) pass=true;
 return pass;
}

//invmass
bool passInvmass(double invmass, double invmass_cuts_low, double invmass_cuts_high)
{
 bool pass=false;
 if(invmass<=invmass_cuts_high && invmass>=invmass_cuts_low) pass=true;
 return pass;
}

//prmary vertex
bool passVertex(double pvz, double pvz_cuts)
{
 bool pass=false;
 if(fabs(pvz)<=pvz_cuts) pass=true;
 return pass;
}

//luminosity
bool passLumosity(double instlum, double instlum_cuts_high, double instlum_cuts_low)
{
 bool pass=false;
 if(instlum>=instlum_cuts_low && instlum<=instlum_cuts_high) pass=true;
 return pass;
}

//efficiency corrections. 3 loops 
double getEffEMIDTrk(double emid_2d, double emid_1d, double trk_2d, double trk_1d)
{
 double eff_rwt;
 if((emid_2d>0.05 && emid_2d<1.5) && (trk_2d>0.05 && trk_2d<1.5))
  eff_rwt = emid_2d*trk_2d;
 else if((emid_2d<=0.05 || emid_2d>=1.5) && (trk_2d>0.05 && trk_2d<1.5))
  eff_rwt = emid_1d*trk_2d;
 else if((emid_2d>0.05 && emid_2d<1.5) && (trk_2d<=0.05 || trk_2d>=1.5))
  eff_rwt = emid_2d*trk_1d;
 else if((emid_2d<=0.05 || emid_2d>=1.5) && (trk_2d<=0.05 || trk_2d>=1.5))
  eff_rwt = emid_1d*trk_1d;

 return eff_rwt;
}

double getEffEMID(double emid_2d, double emid_1d)
{
 double eff_rwt;
 if(emid_2d>0.05 && emid_2d<1.5)
  eff_rwt = emid_2d;
 else 
  eff_rwt = emid_1d;

 return eff_rwt;
}

//energy corrections
//constents

//these functions return the scale factors to the electron 4-momentum
//luminosity correction
double getLumCorr(double deteta, double luminosity, int type)
{
 double global_scale;
  global_scale = 1.0;
 if(type==1001) //iib1 data
  {
   if(fabs(deteta)<=1.1) global_scale = 1.0/(1.00133 - 0.000522655*luminosity);
   else global_scale = 1.0;
  }
 else if(type==1002) //iib2 data
  {
   if(fabs(deteta)<=1.1) global_scale = 1.0/(1.00646 - 23.5438 * exp(0.144927*luminosity-8.58824));
   else if(fabs(deteta)>=1.5) global_scale = pow((1.0/(0.997047 + 0.0116717*exp(-1.0*luminosity))), 2);
   else global_scale = 1.0;
  }
 if(type==1003) //iib3 data
  {
   if(fabs(deteta)<=1.1) global_scale = 1.0/(1.01177 - 0.00172997*luminosity);
   else global_scale = 1.0;
  }
 else if(type==1004) //iib4 data
  {
   if(fabs(deteta)<=1.1)
    {
     if(luminosity<=4.5)
      global_scale = 1.0/(1.01011 - 0.000208855*exp(0.734617*luminosity - 0.00478183));
     else if(luminosity>4.5)
      global_scale = 1.0/(1.00527 - 0.000118274*exp(0.509277*luminosity + 0.0000000678946));
    }
//   else if(fabs(deteta)>=1.5) global_scale = pow((1.0/(0.997047 + 0.0116717*exp(-1.0*luminosity))), 2);
   else global_scale = 1.0;
  }
 else if(type==2) // iib1 mc
  {
   if(fabs(deteta)<=1.1) global_scale = 1.0/(0.993154 + 0.000811637*luminosity);
   else if(fabs(deteta)>=1.5) global_scale = pow((1.0/(1.00603 - 25.0083*exp(-0.184445*luminosity - 8.64968))), 2);
   else global_scale = 1.0;
  }
 else global_scale=1.0;
 return global_scale;
}

//smear parameters
 double mc_fullsmear_ec[14] = {1,1,1,1,1,1,1,0.0101068,0.0150656,0.0100048,0.046914,0.0496904,0.0115974,0.0100203};
 double mc_smear_0003[4] = {0.017547, 0.962147, 7.0, 1.0053};
 double mc_smear_0308[4] = {0.0189464, 0.995784, 7.0, 1.00542};
 double mc_smear_0811[4] = {0.0225773, 0.910745, 7.0, 1.01167};
 double mc_smear_1518[2] = {0.998585, 0.0325467};
 double mc_smear_1821[2] = {1.00175, 0.0178888};
 double mc_smear_2125[2] = {1.00751, 0.0191124};
 double mc_smear_2532[2] = {0.993207, 0.0164855};
/*
 double mc_smear_ccphimod_0005[2] = {0.999445, 0.0207703};
 double mc_smear_ccphimod_0509[2] = {1.00058, 0.0281554 };
 double mc_smear_ccphimod_0911[2] = {1.00052, 0.0380292};
 double mc_smear_ccnonphimod_0005[4] = {1.00536, 0.0403757, 0.9431, 7.0};
 double mc_smear_ccnonphimod_0509[4] = {1.01396, 0.0386896, 0.928585, 7.0};
 double mc_smear_ccnonphimod_0911[4] = {1.00681, 0.00972298, 0.90836, 7.0};
*/
//reprocessed MC
 double mc_smear_ccphimod_0005[2] = {1.00069, 0.0230891};
 double mc_smear_ccphimod_0509[2] = {1.00058, 0.0281554};
 double mc_smear_ccphimod_0911[2] = {1.00106, 0.0345751};
 double mc_smear_ccnonphimod_0005[4] = {1.00066, 0.0347911, 0.908796, 7.0};
 double mc_smear_ccnonphimod_0509[4] = {1.00899, 0.0387548, 0.901872, 7.0};
 double mc_smear_ccnonphimod_0911[4] = {1.00703, 0.0443313, 0.867358, 7.0};

double MCsmear(double deteta, double phimod, double elec_E, TRandom3 *rr)
{
 double global_smear;
 if(fabs(deteta)<0.5)
  {
   if(phimod>0.1 && phimod<0.9)
    global_smear = mc_smear_ccphimod_0005[0]*(1.0+mc_smear_ccphimod_0005[1]*rr->Gaus(0.0,1.0));
   else
    global_smear = mc_smear_ccnonphimod_0005[0]*(1.0 + flipCrystalBall(0, mc_smear_ccnonphimod_0005[1], mc_smear_ccnonphimod_0005[2], mc_smear_ccnonphimod_0005[3], rr));
  }
 else if(fabs(deteta)>0.5 && fabs(deteta)<0.9)
  {
   if(phimod>0.1 && phimod<0.9)
    global_smear = mc_smear_ccphimod_0509[0]*(1.0+mc_smear_ccphimod_0509[1]*rr->Gaus(0.0,1.0));
   else
    global_smear = mc_smear_ccnonphimod_0509[0]*(1.0 + flipCrystalBall(0, mc_smear_ccnonphimod_0509[1], mc_smear_ccnonphimod_0509[2], mc_smear_ccnonphimod_0509[3], rr));
  }
 else if(fabs(deteta)>0.9 && fabs(deteta)<1.1)
  {
   if(phimod>0.1 && phimod<0.9)
    global_smear = mc_smear_ccphimod_0911[0]*(1.0+mc_smear_ccphimod_0911[1]*rr->Gaus(0.0,1.0));
   else
    global_smear = mc_smear_ccnonphimod_0911[0]*(1.0 + flipCrystalBall(0, mc_smear_ccnonphimod_0911[1], mc_smear_ccnonphimod_0911[2], mc_smear_ccnonphimod_0911[3], rr));
  }
 else
  global_smear = 1.0;

 return global_smear;
}

double getMCfullSmear(double deteta, double phimod, double elec_E, TRandom3 *rr)
{
  double global_scale;
  if(fabs(deteta)<=1.1)
   {
    if(phimod>0.1 && phimod<0.9)
     {
      if(fabs(deteta)<0.5)
       global_scale = 1.0*(1.0+0.0218499*rr->Gaus(0.0,1.0));
      else if(fabs(deteta)>0.5 && fabs(deteta)<0.9)
       global_scale = 1.0*(1.0+0.02635*rr->Gaus(0.0,1.0));
      else if(fabs(deteta)>0.9 && fabs(deteta)<1.1)
       global_scale = 1.0*(1.0+0.0255879*rr->Gaus(0.0,1.0));
     }
    else if(phimod<=0.1 || phimod>=0.9)
//     global_scale = 1.0156*(1.0 + flipCrystalBall(0, 0.0342758, 0.759312, 7.0, rr));
     global_scale = 1.013*(1.0 + flipCrystalBall(0, 0.0394656, 0.840258, 7.0, rr));
   }
  else if(fabs(deteta)>=1.5)
   {
    if(fabs(deteta)>1.5 && fabs(deteta)<1.7) global_scale = mc_fullsmear_ec[0]*(1.0 + mc_fullsmear_ec[0+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>1.7 && fabs(deteta)<1.9) global_scale = mc_fullsmear_ec[1]*(1.0 + mc_fullsmear_ec[1+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>1.9 && fabs(deteta)<2.1) global_scale = mc_fullsmear_ec[2]*(1.0 + mc_fullsmear_ec[2+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>2.1 && fabs(deteta)<2.4) global_scale = mc_fullsmear_ec[3]*(1.0 + mc_fullsmear_ec[3+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>2.4 && fabs(deteta)<2.6) global_scale = mc_fullsmear_ec[4]*(1.0 + mc_fullsmear_ec[4+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>2.6 && fabs(deteta)<2.8) global_scale = mc_fullsmear_ec[5]*(1.0 + mc_fullsmear_ec[5+7]*rr->Gaus(0.0,1.0));
    if(fabs(deteta)>2.8 && fabs(deteta)<3.2) global_scale = mc_fullsmear_ec[6]*(1.0 + mc_fullsmear_ec[6+7]*rr->Gaus(0.0,1.0));

   }
  else  global_scale = 1;
    return global_scale;
}

double getMCsimpleSmear(double deteta, double phimod, double elec_E, TRandom3 *rr)
{
  double global_scale;
  if(fabs(deteta)<=1.1)
   {
    if(phimod>0.1 && phimod<0.9)
      global_scale = 1.0*(1.0+0.0244849*rr->Gaus(0.0,1.0));
//      global_scale = 0.999817*(1.0+0.0244849*rr->Gaus(0.0,1.0));
    else if(phimod<=0.1 || phimod>=0.9)
     global_scale = 1.00878*(1.0 + flipCrystalBall(0, 0.0403415, 0.897159, 7.0, rr)); 
   }
  else if(fabs(deteta)>=1.5)
   global_scale = 1.0*(1.0+0.0253267*rr->Gaus(0.0,1.0));
//   global_scale = 0.998556*(1.0+0.0247992*rr->Gaus(0.0,1.0));

  else  global_scale = 1;
    return global_scale;
}

double getMCSpecialSmear(double deteta, double phimod, double elec_E, TRandom3 *rr)
{
  double global_scale;
  double etasmear;
  if(fabs(deteta)<=1.1)
   {
    if(phimod>0.1 && phimod<0.9)
      global_scale = 1.0*(1.0+0.0244849*rr->Gaus(0.0,1.0));
    else if(phimod<=0.1 || phimod>=0.9)
     global_scale = 1.00878*(1.0 + flipCrystalBall(0, 0.0403415, 0.897159, 7.0, rr));
   }
  else if(fabs(deteta)>=1.5)
   {
    etasmear = (fabs(deteta)-1.5) / 1.7;
    global_scale = 1.0*(1.0+(0.0253267 + etasmear*0.05)*rr->Gaus(0.0,1.0));
   }
  else  global_scale = 1;
    return global_scale;
}

double getCBSmearCorr(double deteta, double phimod, double elec_E, TRandom3 *rr)
{
 double global_scale;
 if(fabs(deteta)<=1.1)//in CC
  {
   if(fabs(deteta)<0.3)//0.0~0.3
     global_scale = mc_smear_0003[3]*(1.0+flipCrystalBall(0, mc_smear_0003[0], mc_smear_0003[1], mc_smear_0003[2], rr));
   else if(fabs(deteta)>0.3 && fabs(deteta)<0.8)//0.3~0.8
     global_scale = mc_smear_0308[3]*(1.0+flipCrystalBall(0, mc_smear_0308[0], mc_smear_0308[1], mc_smear_0308[2], rr));
   else if(fabs(deteta)>0.8)//0.8~1.1
     global_scale = mc_smear_0811[3]*(1.0+flipCrystalBall(0, mc_smear_0811[0], mc_smear_0811[1], mc_smear_0811[2], rr));
  }
 else if(fabs(deteta)>=1.5)//in EC 
  {
   if(fabs(deteta)<1.8)//1.5~1.8
    global_scale = mc_smear_1518[0]*(1.0+mc_smear_1518[1]*rr->Gaus(0.0,1.0));
   else if(fabs(deteta)>1.8 && fabs(deteta)<2.1)//1.8~2.1
    global_scale = mc_smear_1821[0]*(1.0+mc_smear_1821[1]*rr->Gaus(0.0,1.0));
   else if(fabs(deteta)>2.1 && fabs(deteta)<2.5)//2.1~2.5
    global_scale = mc_smear_2125[0]*(1.0+mc_smear_2125[1]*rr->Gaus(0.0,1.0));
   else if(fabs(deteta)>2.5)
    global_scale = mc_smear_2532[0]*(1.0+mc_smear_2532[1]*rr->Gaus(0.0,1.0));
  }
 else
  global_scale = 1.0;

 return global_scale;
}

//eta correction factors
  double cccc_scale_data[10];
  double cccc_scale_data_unc[10];
  double data_eta_corr[14];
  double data_eta_corr_sc[14];
  double pp_data_1[24];
  double pp_data_0[24];

//fot reprocessed iib2 data
  double cccc_scale_data_iib2[10] = {0.999066, 0.998968, 0.998931, 1.00217, 1.00143, 1.00026, 0.999768, 0.998966, 0.9981, 0.999649};
  double cccc_scale_data_unc_iib2[10] = {0.0016,0.0016,-0.00399999,0.00999983,-0.00399999,0.389418,0.00999983,0.00999983,0,0};
  double data_eta_corr_iib2[14] = {1.099, 1.04248, 1.03862, 1.04397, 1.00724, 1.02782, 0.999865, 1.02289, 1.00139, 1.01159, 1.02412, 1.05592, 1.06823, 1.09992};
  double data_eta_corr_sc_iib2[14] = {0.35866, 0.0320008, 0.0171978, 0.0126273, 0.185416, 0.0120191, 0.00818091, 0.00437863, 0.00540245, 0.0137124, -0.0365292, 0.0230209, 0.0896632, 0.0117967};
  double pp_data_0_iib2[24] = {281.061, 227.547, 194.44, 174.915, 161.904, 147.785, 136.374, 127.56, 116.13, 107.912, 99.807, 91.8785, 89.6821, 98.2954, 105.47, 114.181, 124.273, 135.355, 144.722, 156.247, 171.17, 189.566, 221.85, 280.801};
  double pp_data_1_iib2[24] = {-291.638, -232.658, -193.929, -174.151, -158.247, -146.957, -136.674, -126.487, -114.556, -107.35, -99.6225, -91.4833, -90.4498, -99.342, -105.928, -114.205, -123.797, -135.684, -144.44, -155.836, -173.82, -189.889, -225.968, -288.994};
   
//for reprocessed iib1 data
  double cccc_scale_data_iib1[10] = {1.00047, 0.998044, 1.00236, 0.999975, 0.999018, 1.00148, 1.00069, 0.998503, 0.999035, 0.998608};
  double cccc_scale_data_unc_iib1[10] = {-0.00399999,0.00639996,0.389418,0.129634,0.00999983,0.00999983,0.0399893,0.0399893,0.0399893,0};
  double data_eta_corr_iib1[14] = {1.07128, 1.02097, 1.08639, 1.00841, 1.01674, 1.00654, 1.0174, 0.989929, 1.00703, 1.01643, 1.01946, 1.05807, 1.05292, 1.03907};
  double data_eta_corr_sc_iib1[14] = {0.0105628, -0.000935459, 0.00122773, -0.00996432, -0.0724938, -0.0186422, -0.0113631, -0.0109964, -0.0122537, -0.0538925, -0.0424924, -0.00102821, 8.83937e-05, 0.209675};
  double pp_data_0_iib1[24] = {299.892, 218.013, 192.446, 173.391, 160.355, 147.688, 137.083, 124.321, 117.257, 108.537, 99.3014, 91.8642, 90.1131, 97.4832, 106.208, 114.295, 127.536, 133.968, 145.851, 155.552, 166.941, 194.825, 223.23, 276.238};
  double pp_data_1_iib1[24] = {-309.908, -222.75, -194.183, -174.51, -157.152, -147.353, -138.271, -124.456, -116.549, -108.435, -99.6304, -91.9857, -90.4123, -98.0615, -106.284, -113.916, -127.462, -133.408, -144.939, -156.447, -169.725, -194.085, -226.111, -280.675};

//fot reprocessed iib3 data
  double cccc_scale_data_iib3[10] = {0.998306, 0.99763, 0.999582, 0.99839, 0.999228, 0.995889, 0.998991, 1.00441, 1.00307, 1.00638};
  double cccc_scale_data_unc_iib3[10] = {0.00952034, -0.041855, 0.00312895, 0.00340418, -0.00943444, 0.00525045, 0.00105568, 0.00570976, 0.00155138, 0.00109903};
  double data_eta_corr_iib3[14] = {1.08671, 1.04629, 1.04701, 1.00694, 1.01517, 1.00381, 1.00192, 1.00858, 0.99555, 1.00374, 1.01074, 1.01363, 1.08658, 1.09888};
  double data_eta_corr_sc_iib3[14] = {0.0433933, 0.0157058, -0.174665, 0.112919, 0.0188607, 0.0138292, 0.0110359, 0.0104814, -0.00900813, -0.0905204, 0.00952456, 0.00212531, -0.0115311, 0.00632724};
  double pp_data_0_iib3[24] = {291.8, 233.356, 196.565, 173.968, 159.906, 146.045, 138.448, 128.441, 117.692, 109.114, 101.036, 93.0447, 90.3534, 98.8543, 104.573, 114.001, 122.749, 135.479, 144.89, 156.873, 173.73, 192.741, 215.75, 284.08};
  double pp_data_1_iib3[24] = {-293.593, -233.967, -198.041, -175.91, -161.072, -146.799, -138.879, -128.012, -118.84, -109.46, -101.413, -92.6937, -90.1836, -99.1823, -104.255, -114.999, -123.385, -136.947, -145.064, -156.988, -175.215, -193.321, -217.719, -286.238};

//for reprocessed iib4 data
  double cccc_scale_data_iib4[10] = {0.993254, 0.996891, 1.00177, 0.995281, 1.00135, 0.999949, 0.998254, 1.00114, 1.00626, 1.00205};
  double cccc_scale_data_unc_iib4[10] = {-0.00376999, 0.00821741, -0.00346399, -0.002981, 0.00624496, -0.002015, -0.001532, -0.001049, -0.000566, 0.0002075};
  double data_eta_corr_iib4[14] = {1.07769, 1.03647, 1.03156, 1.01109, 1.00679, 1.00181, 1.00428, 1.02021, 1.0016, 1.01763, 1.01994, 1.02235, 1.05846, 1.09598};
  double data_eta_corr_sc_iib4[14] = {0.230958, 0.0185087, 0.0247558, 0.0364307, 0.0279797, 0.016598, 0.0364777, 0.00239776, 0.0019052, -0.135324, -0.00933142, 0.00851286, 0.00399823, 0.131036};
  double pp_data_0_iib4[24] = {271.49, 224.598, 192.751, 175.194, 157.533, 150.145, 138.001, 128.416, 117.292, 109.344, 99.9247, 91.9207, 90.1894, 97.985, 106.408, 116.284, 123.804, 136.76, 145.845, 157.396, 171.6, 191.202, 218.544, 275.713};
  double pp_data_1_iib4[24] = {-275.832, -227.068, -194.117, -175.956, -159.037, -151.15, -138.974, -128.822, -117.676, -109.996, -100.298, -92.1377, -90.4059, -98.35, -106.459, -117.143, -124.623, -137.466, -147.018, -158.951, -173.756, -192.916, -220.03, -279.256};

//for reprocessed iib1 MC, full measurement
  double cccc_scale_mc[10] = {0.994762, 0.999207, 0.99981, 1.00153, 1.00152, 1.00185, 1.00099, 0.999894, 0.99862, 0.994357};
  double cccc_scale_mc_unc[10] = {0,0,0,0,0,0,0,0,0,0};
  double mc_eta_corr[14] = {1.11548, 1.06042, 1.03654, 1.0329, 1.02539, 1.02883, 1.02985, 1.02713, 1.02444, 1.0208, 1.02678, 1.03834, 1.06583, 1.10669};
  double mc_eta_corr_sc[14] = {0.0136402, -8.33506e-05, 0.00637541, 0.00281486, 0.00654599, 0.00734244, 0.0160203, 0.0205874, 0.0245514, 0.0373316, 0.000473161, 0.00817936, -0.000279666, -0.110346};
//  double mc_eta_corr[14] = {1.11571, 1.06178, 1.03887, 1.03317, 1.02571, 1.02849, 1.03002, 1.02715, 1.0245, 1.02028, 1.02776, 1.04029, 1.06636, 1.10676};
//  double mc_eta_corr_sc[14] = {0.0142861, 0.000142033, 0.00643842, 0.00350262, 0.00704256, 0.00740386, 0.00760068, 0.00387024, 0.011181, 0.00319919, 7.96512e-05, 0.00723679, 0.000705815, 0.000314528};
  double pp_mc_0[24] = {277.136, 225.052, 191.311, 169.812, 157.728, 145.676, 134.891, 124.83, 115.207, 106.475, 98.7842, 91.1501, 90.0662, 98.3232, 105.074, 113.807, 122.744, 133.016, 143.722, 155.442, 165.811, 188.05, 220.5, 270.733};
  double pp_mc_1[24] = {-285.488, -229.194, -193.335, -170.637, -158.245, -145.798, -134.653, -124.296, -114.638, -105.91, -98.2064, -90.9118, -90.0292, -97.7276, -104.53, -113.227, -122.828, -132.695, -143.77, -155.798, -166.435, -190.147, -224.434, -278.448};

double getEtaCorr(double deteta, double elecE, double lumi, int type)
{
 double global_scale;
 int parRegion;
 int EtaRegion;
 EtaRegion = (int)((deteta+3.2)/0.1);
 if(type==2) //for MC
  {// MC
   if(fabs(deteta)<=1.1) //for CC electron
    {
     if(EtaRegion>=21 && EtaRegion<23) global_scale = (cccc_scale_mc[0]*elecE + cccc_scale_mc_unc[0])/elecE;
     if(EtaRegion>=23 && EtaRegion<25) global_scale = (cccc_scale_mc[1]*elecE + cccc_scale_mc_unc[1])/elecE;
     if(EtaRegion>=25 && EtaRegion<27) global_scale = (cccc_scale_mc[2]*elecE + cccc_scale_mc_unc[2])/elecE;
     if(EtaRegion>=27 && EtaRegion<29) global_scale = (cccc_scale_mc[3]*elecE + cccc_scale_mc_unc[3])/elecE;
     if(EtaRegion>=29 && EtaRegion<32) global_scale = (cccc_scale_mc[4]*elecE + cccc_scale_mc_unc[4])/elecE;
     if(EtaRegion>=32 && EtaRegion<35) global_scale = (cccc_scale_mc[5]*elecE + cccc_scale_mc_unc[5])/elecE;
     if(EtaRegion>=35 && EtaRegion<37) global_scale = (cccc_scale_mc[6]*elecE + cccc_scale_mc_unc[6])/elecE;
     if(EtaRegion>=37 && EtaRegion<39) global_scale = (cccc_scale_mc[7]*elecE + cccc_scale_mc_unc[7])/elecE;
     if(EtaRegion>=39 && EtaRegion<41) global_scale = (cccc_scale_mc[8]*elecE + cccc_scale_mc_unc[8])/elecE;
     if(EtaRegion>=41 && EtaRegion<43) global_scale = (cccc_scale_mc[9]*elecE + cccc_scale_mc_unc[9])/elecE;
    }
   else if(fabs(deteta)>=1.5)
    {

     if(EtaRegion<4) parRegion=0;
     if(EtaRegion>=4 && EtaRegion<6) parRegion=1;
     if(EtaRegion>=6 && EtaRegion<8) parRegion=2;
     if(EtaRegion==8 || EtaRegion==9 || EtaRegion==10) parRegion=3;
     if(EtaRegion==11 || EtaRegion==12) parRegion=4;
     if(EtaRegion==13 || EtaRegion==14) parRegion=5;
     if(EtaRegion==15 || EtaRegion==16) parRegion=6;
     if(EtaRegion==47 || EtaRegion==48) parRegion=7;
     if(EtaRegion==49 || EtaRegion==50) parRegion=8;
     if(EtaRegion==51 || EtaRegion==52) parRegion=9;
     if(EtaRegion==53 || EtaRegion==54 || EtaRegion==55) parRegion=10;
     if(EtaRegion>=56 && EtaRegion<58) parRegion=11;
     if(EtaRegion>=58 && EtaRegion<60) parRegion=12;
     if(EtaRegion>=60 && EtaRegion<64) parRegion=13;

     if( EtaRegion<4 ) //-3.2 ~ -2.8
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[0]*mc_eta_corr[parRegion]+pp_mc_0[0]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=4 && EtaRegion<6) //-2.8 ~ -2.6
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[1]*mc_eta_corr[parRegion]+pp_mc_0[1]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=6 && EtaRegion<8) //-2.6 ~ -2.4
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[2]*mc_eta_corr[parRegion]+pp_mc_0[2]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=8 && EtaRegion<17) //-2.4 ~ -1.5
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[EtaRegion-5]*mc_eta_corr[parRegion]+pp_mc_0[EtaRegion-5]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=47 && EtaRegion<56) //+1.5 ~ +2.4
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[EtaRegion-35]*mc_eta_corr[parRegion]+pp_mc_0[EtaRegion-35]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=56 && EtaRegion<58) //+2.4 ~ +2.6
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[21]*mc_eta_corr[parRegion]+pp_mc_0[21]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=58 && EtaRegion<60) //+2.6 ~ +2.8
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[22]*mc_eta_corr[parRegion]+pp_mc_0[22]+mc_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=60 && EtaRegion<64) //+2.8 ~ +3.2
       global_scale = (mc_eta_corr[parRegion]*elecE+pp_mc_1[23]*mc_eta_corr[parRegion]+pp_mc_0[23]+mc_eta_corr_sc[parRegion])/elecE;

//     global_scale = 1.0;
    }
  }// MC
 else if(type==1001 || type==1002 || type==1003 || type==1004) //for data
  {// data
   if(type==1001)
    {
      for(int ipart=0;ipart<10;ipart++)
       {
        cccc_scale_data[ipart] = cccc_scale_data_iib1[ipart];
//        cccc_scale_data_unc[ipart] = cccc_scale_data_unc_iib1[ipart];
        cccc_scale_data_unc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<14;ipart++)
       {
        data_eta_corr[ipart] = data_eta_corr_iib1[ipart];
//        data_eta_corr_sc[ipart] = data_eta_corr_sc_iib1[ipart];
        data_eta_corr_sc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<24;ipart++)
       {
        pp_data_1[ipart] = pp_data_1_iib1[ipart];
        pp_data_0[ipart] = pp_data_0_iib1[ipart];
       }
    }
   else if(type==1002)
    {
      for(int ipart=0;ipart<10;ipart++)
       {
        cccc_scale_data[ipart] = cccc_scale_data_iib2[ipart];
//        cccc_scale_data_unc[ipart] = cccc_scale_data_unc_iib2[ipart];
        cccc_scale_data_unc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<14;ipart++)
       {
        data_eta_corr[ipart] = data_eta_corr_iib2[ipart];
        data_eta_corr_sc[ipart] = data_eta_corr_sc_iib2[ipart];
//        data_eta_corr_sc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<24;ipart++)
       {
        pp_data_1[ipart] = pp_data_1_iib2[ipart];
        pp_data_0[ipart] = pp_data_0_iib2[ipart];
       }
    }
   else if(type==1003)
    {
      for(int ipart=0;ipart<10;ipart++)
       {
        cccc_scale_data[ipart] = cccc_scale_data_iib3[ipart];
        cccc_scale_data_unc[ipart] = cccc_scale_data_unc_iib3[ipart];
//        cccc_scale_data_unc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<14;ipart++)
       {
        data_eta_corr[ipart] = data_eta_corr_iib3[ipart];
        data_eta_corr_sc[ipart] = data_eta_corr_sc_iib3[ipart];
//        data_eta_corr_sc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<24;ipart++)
       {
        pp_data_1[ipart] = pp_data_1_iib3[ipart];
        pp_data_0[ipart] = pp_data_0_iib3[ipart];
       }
    }
   else if(type==1004)
    {
      for(int ipart=0;ipart<10;ipart++)
       {
        cccc_scale_data[ipart] = cccc_scale_data_iib4[ipart];
//        cccc_scale_data_unc[ipart] = cccc_scale_data_unc_iib4[ipart];
        cccc_scale_data_unc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<14;ipart++)
       {
        data_eta_corr[ipart] = data_eta_corr_iib4[ipart];
        data_eta_corr_sc[ipart] = data_eta_corr_sc_iib4[ipart];
//        data_eta_corr_sc[ipart] = 0.0;
       }
      for(int ipart=0;ipart<24;ipart++)
       {
        pp_data_1[ipart] = pp_data_1_iib4[ipart];
        pp_data_0[ipart] = pp_data_0_iib4[ipart];
       } 
    } 
   else cout<<"what the ... ..."<<endl;

   if(fabs(deteta)<=1.1) //for CC electron
    {
     if(EtaRegion>=21 && EtaRegion<23) global_scale = (cccc_scale_data[0]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=23 && EtaRegion<25) global_scale = (cccc_scale_data[1]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=25 && EtaRegion<27) global_scale = (cccc_scale_data[2]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=27 && EtaRegion<29) global_scale = (cccc_scale_data[3]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=29 && EtaRegion<32) global_scale = (cccc_scale_data[4]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=32 && EtaRegion<35) global_scale = (cccc_scale_data[5]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=35 && EtaRegion<37) global_scale = (cccc_scale_data[6]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=37 && EtaRegion<39) global_scale = (cccc_scale_data[7]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=39 && EtaRegion<41) global_scale = (cccc_scale_data[8]*elecE + cccc_scale_data_unc[0])/elecE;
     if(EtaRegion>=41 && EtaRegion<43) global_scale = (cccc_scale_data[9]*elecE + cccc_scale_data_unc[0])/elecE;
    }
   else if(fabs(deteta)>=1.5)
    {
//     global_scale = 1.0;

     if(EtaRegion<4) parRegion=0;
     if(EtaRegion>=4 && EtaRegion<6) parRegion=1;
     if(EtaRegion>=6 && EtaRegion<8) parRegion=2;
     if(EtaRegion==8 || EtaRegion==9 || EtaRegion==10) parRegion=3;
     if(EtaRegion==11 || EtaRegion==12) parRegion=4;
     if(EtaRegion==13 || EtaRegion==14) parRegion=5;
     if(EtaRegion==15 || EtaRegion==16) parRegion=6;
     if(EtaRegion==47 || EtaRegion==48) parRegion=7;
     if(EtaRegion==49 || EtaRegion==50) parRegion=8;
     if(EtaRegion==51 || EtaRegion==52) parRegion=9;
     if(EtaRegion==53 || EtaRegion==54 || EtaRegion==55) parRegion=10;
     if(EtaRegion>=56 && EtaRegion<58) parRegion=11;
     if(EtaRegion>=58 && EtaRegion<60) parRegion=12;
     if(EtaRegion>=60 && EtaRegion<64) parRegion=13;

     if( EtaRegion<4 ) //-3.2 ~ -2.8
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[0]*data_eta_corr[parRegion]+pp_data_0[0]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=4 && EtaRegion<6) //-2.8 ~ -2.6
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[1]*data_eta_corr[parRegion]+pp_data_0[1]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=6 && EtaRegion<8) //-2.6 ~ -2.4
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[2]*data_eta_corr[parRegion]+pp_data_0[2]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=8 && EtaRegion<17) //-2.4 ~ -1.5
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[EtaRegion-5]*data_eta_corr[parRegion]+pp_data_0[EtaRegion-5]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=47 && EtaRegion<56) //+1.5 ~ +2.4
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[EtaRegion-35]*data_eta_corr[parRegion]+pp_data_0[EtaRegion-35]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=56 && EtaRegion<58) //+2.4 ~ +2.6
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[21]*data_eta_corr[parRegion]+pp_data_0[21]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=58 && EtaRegion<60) //+2.6 ~ +2.8
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[22]*data_eta_corr[parRegion]+pp_data_0[22]+data_eta_corr_sc[parRegion])/elecE;
     else if( EtaRegion>=60 && EtaRegion<64) //+2.8 ~ +3.2
       global_scale = (data_eta_corr[parRegion]*elecE+pp_data_1[23]*data_eta_corr[parRegion]+pp_data_0[23]+data_eta_corr_sc[parRegion])/elecE;

    }
   else
    global_scale = 1.0;
  }// data
 else global_scale = 1.0;
 return global_scale;
}


//QCD shape correction
double getQCDshapeCorr(double pt, double deteta)
{
 double rwt;
 double jet_rwt;
 double elec_rwt;

 rwt = 1.0;
 if(fabs(deteta)<1.1)
  {
   jet_rwt = 0.667118 - 0.0893544*exp(-0.0875087*(pt - 31.3103));
//   elec_rwt = 0.0353914*exp(-0.135335*(pt-32.2499));
   elec_rwt = 0.0353914*exp(-0.135335*(pt-32.2499)) + 0.0479033*exp(0.00218316*pt+1.79851);
  }
 else if(fabs(deteta)>1.5)
  {
   jet_rwt = 0.548977 - 0.0557609*exp(-0.0810507*(pt - 34.4091));
//   elec_rwt = 0.0155555*exp(-0.173519*(pt - 33.1874));
   elec_rwt = 0.0155555*exp(-0.173519*(pt - 33.1874)) + 0.0420803*exp(0.0013439*pt+2.08121);
  }

 rwt = elec_rwt / jet_rwt;
// if(fabs(deteta)<1.1) cout<<pt<<"  "<<elec_rwt<<"  "<<jet_rwt<<"  "<<rwt<<endl;
 return rwt;
}

//re-calculate angle
double recPolarAngle(double invtxz, double inrevtxz, double deteta, double phyeta)
{
  double theta_D, theta_E, theta_R;
  double tanth_D, tanth_E, tanth_R;
  double OB, AB;

  tanth_D = atan(2*exp(-1.0*deteta)/(1 - exp(-2.0*deteta)));
  tanth_E = atan(2*exp(-1.0*phyeta)/(1 - exp(-2.0*phyeta)));
  OB = invtxz / (1-(tanth_D/tanth_E));
  AB = tanth_D * OB;
  tanth_R = AB / (OB - inrevtxz);
  theta_R = atan(tanth_R);
  return theta_R;
}

//cut the edge of every eta bin
bool passCutEtaEdge(double deteta)
{
 bool finalpass=true;
 double edge_eta[49] = {-3.2,-2.8,-2.6,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.6,2.8,3.2};
 for(int ibin=0;ibin<49;ibin++)
  {
   if(fabs(deteta-edge_eta[ibin])<0.01) finalpass=false;
  }

 return finalpass;
}

//re-calculate phyeta
double recPhyEta(double invtxz, double inrevtxz, double deteta, double phyeta)
{
  double theta_D, theta_E, theta_R;
  double tanth_D, tanth_E, tanth_R;
  double OB, AB;
  double reEta;

  tanth_D = atan(2*exp(-1.0*deteta)/(1 - exp(-2.0*deteta)));
  tanth_E = atan(2*exp(-1.0*phyeta)/(1 - exp(-2.0*phyeta)));
  OB = invtxz / (1-(tanth_D/tanth_E));
  AB = tanth_D * OB;
  tanth_R = AB / (OB - inrevtxz);

  theta_R = atan(tanth_R);
  reEta = -1.0 * log(theta_R/0.5);
  return reEta;
}
