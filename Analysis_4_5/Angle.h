#include <iostream>
#include <fstream>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1D.h"
 
  double func_cos_theta( TLorentzVector elec0, TLorentzVector elec1, double Q, double Qt)
{
  
  double cos_theta =   1  /(Q*sqrt(Q*Q+Qt*Qt))*  ((elec0.E()+elec0.Pz())*(elec1.E()-elec1.Pz())-
                 (elec0.E()-elec0.Pz())*(elec1.E()+elec1.Pz()));
  return cos_theta;
} 



  double func_phi( TLorentzVector elec0, TLorentzVector elec1)
{
  TLorentzVector dilepton = elec0 + elec1 ;
  TVector3 boost          = -dilepton.BoostVector();

  TLorentzVector el(elec0.Px(), elec0.Py(), elec0.Pz(), elec0.E());
  TLorentzVector po(elec1.Px(), elec1.Py(), elec1.Pz(), elec1.E());
  el.Boost(boost);
  po.Boost(boost);

  TLorentzVector lab_z1(0.,0., 980.,980.);
  TLorentzVector lab_z2(0.,0.,-980.,980.);
  lab_z1.Boost(boost);
  lab_z2.Boost(boost);

  TVector3 boson_frame  = (lab_z1.Vect().Unit()-lab_z2.Vect().Unit()).Unit();
  TVector3 yplant_frame = ((lab_z1.Vect().Unit()).Cross(lab_z2.Vect().Unit())).Unit();
  TVector3 xplant_frame = (yplant_frame.Cross(boson_frame)).Unit();
  TVector3 elec_vec     = el.Vect();
  double phi = atan2(elec_vec.Dot(yplant_frame), elec_vec.Dot(xplant_frame)); 
  return phi;
}
 
