/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include <map>
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
//  private:
 protected:
  std::shared_ptr<Branches6> _data;

  double _beam_energy = E1D_E0;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;

  std::unique_ptr<TLorentzVector> _boosted_gamma;
  std::unique_ptr<TLorentzVector> _boosted_prot;
  std::unique_ptr<TLorentzVector> _boosted_pip;
  std::unique_ptr<TLorentzVector> _boosted_pim;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  float _MM = std::nanf("-99");
  float _MM2 = std::nanf("-99");

  float _W = std::nanf("-99");
  float _Q2 = std::nanf("-99");

  float _inv_Ppip = NAN;
  float _inv_Ppim = NAN;
  float _inv_pip_pim = NAN;
  float _W_P2pi = NAN;

  void SetElec();

 public:
  Reaction(const std::shared_ptr<Branches6> data);
  ~Reaction();

  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  void CalcMissMass();
  float MM();
  float MM2();
  
  float _phi_gamma = NAN;
  float _phi_prot = NAN;
  float _phi_pip = NAN;
  float _phi_pim = NAN;

  float _alpha_ppip_pipim = NAN;
  float _alpha_pippim_pipf = NAN;
  float _alpha_ppim_pipip = NAN;

  void boost();
  bool _is_boosted = false;

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  inline float weight() {
         return 1.0;
  }

  float inv_Ppip();
  float inv_Ppim();
  float inv_pip_pim();
  float w_P2pi_rec();

  void W_2pi_P();
  void invMassPpip();
  void invMassPpim();
  void invMasspippim();

  float prot_theta();
  float pip_theta();
  float pim_theta();

  void AlphaCalc();

  float gamma_Phi();
  float prot_Phi();
  float pip_Phi();
  float pim_Phi();
  float prot_Phi_lab();
  float pip_Phi_lab();
  float pim_Phi_lab();
  float pim_Phi_lab_measured();

  float alpha_ppip_pipim();
  float alpha_pippim_pipf();
  float alpha_ppim_pipip();

  inline bool MM_cut() {
    return  // abs(Reaction::MM2()) < 0.03 // MM2_cut;
        (Reaction::MM2() < 0.08 && Reaction::MM2() > -0.06);
  }

  inline bool TwoPion() {
    return ((_numPip == 1 && _numPim == 1) && (_hasE && !_hasP && _hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool ProtonPim() {
    return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && !_hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool TwoPion_missingPim() {
    return ((_numPip == 1 && _numProt == 1) && (_hasE && _hasPip && _hasP ));
  }
  inline bool SingleP() {
    return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool NeutronPip() {
    return ((_numPip == 1 && _numNeutral == 1) &&
            (_hasE && !_hasP && _hasPip && !_hasPim && _hasNeutron && !_hasOther));
  }

  inline TLorentzVector e_mu() { return *_beam; }
  inline TLorentzVector e_mu_prime() { return *_elec; }
  inline TLorentzVector gamma() { return *_gamma; }
};

#endif
