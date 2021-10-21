/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(std::shared_ptr<Branches6> data) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  if (getenv("BEAM_E") != NULL) _beam_energy = atof(getenv("BEAM_E"));
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();
  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->id(i) == NEUTRON)
    SetNeutron(i);
  else {
    _numOther++;
    _hasOther = true;
    _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass_map[_data->id(i)]);
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  *mm += (*_beam - *_elec + *_target);
  if (TwoPion_missingPim()) {
    *mm -= *_pip;
    *mm -= *_prot;

    _MM = mm->M();
    _MM2 = mm->M2();

  } else if (TwoPion()) {
    *mm -= *_prot;
    *mm -= *_pip;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (ProtonPim()) {
    *mm -= *_prot;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (SingleP()) {
    *mm -= *_prot;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}

void Reaction::invMassPpim() {
  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_prot;
  *inv_Ppim += (*_gamma + *_target - *_prot - *_pip);
  if (TwoPion_missingPim()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_pip;
  *inv_pip_pim += (*_gamma + *_target - *_prot - *_pip);
  if (TwoPion_missingPim()) _inv_pip_pim = inv_pip_pim->M();
}
void Reaction::invMassPpip() {
  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_prot;
  *inv_Ppip += *_pip;
  if (TwoPion_missingPim()) _inv_Ppip = inv_Ppip->M();
}

void Reaction::W_2pi_P() {
  auto W_P2pi = std::make_unique<TLorentzVector>();
  *W_P2pi += *_prot;
  *W_P2pi += *_pip;
  *W_P2pi += (*_gamma + *_target - *_prot - *_pip);

  if (TwoPion_missingPim()) _W_P2pi = W_P2pi->M();

}  
float Reaction::inv_Ppip() {
  if (_inv_Ppip != _inv_Ppip) invMassPpip();
  return _inv_Ppip;
}
float Reaction::inv_Ppim() {
  if (_inv_Ppim != _inv_Ppim) invMassPpim();
  return _inv_Ppim;
}
float Reaction::inv_pip_pim() {
  if (_inv_pip_pim != _inv_pip_pim) invMasspippim();
  return _inv_pip_pim;
}
float Reaction::w_P2pi_rec() {
  if (_W_P2pi != _W_P2pi) W_2pi_P();
  return _W_P2pi;
}

void Reaction::boost() {
  _is_boosted = true;
  _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
  _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
  _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  // _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  TRotation rot;
  _boosted_gamma->Transform(rot);
  float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));
  TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
  ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot
  _boosted_prot->Transform(rot);
  _boosted_prot->Boost(0, 0, -beta_1);
  _boosted_pip->Transform(rot);
  _boosted_pip->Boost(0, 0, -beta_1);
  _boosted_pim->Transform(rot);
  _boosted_pim->Boost(0, 0, -beta_1);
  _boosted_gamma->Boost(0, 0, -beta_1);

  // _boosted_pim_measured->Transform(rot);
  // _boosted_pim_measured->Boost(0, 0, -beta_1);
  // -beta ko value (0.5 to -0.5 huda
  // samma value aauchha nattra aaudyna)

}

float Reaction::prot_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_prot->Theta() * 180.0 / PI;
  // else
  return NAN;
}
float Reaction::pip_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_pip->Theta() * 180.0 / PI;
  // else
  return NAN;
}
float Reaction::pim_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_pim->Theta() * 180.0 / PI;
  // else
  return NAN;
}

float Reaction::gamma_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_gamma->Phi() > 0)
      return _boosted_gamma->Phi() * 180 / PI;
    else if (_boosted_gamma->Phi() < 0)
      return (_boosted_gamma->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::prot_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_prot->Phi() > 0)
      return _boosted_prot->Phi() * 180 / PI;
    else if (_boosted_prot->Phi() < 0)
      return (_boosted_prot->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::pip_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_pip->Phi() > 0)
      return _boosted_pip->Phi() * 180 / PI;
    else if (_boosted_pip->Phi() < 0)
      return (_boosted_pip->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pim_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_pim->Phi() > 0)
      return _boosted_pim->Phi() * 180 / PI;
    else if (_boosted_pim->Phi() < 0)
      return (_boosted_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

void Reaction::AlphaCalc() {
  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm;  // proton initial pim
  float alpha_PIpPIm_pipf;
  float alpha_PPIm_piPIp;

  if (!_is_boosted) boost();

  // 1 this one is used for α[π−]
  a_gamma = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * V3_anti_z),
                               2)));  // V3_anti_z(0,0,-1);
  b_gamma = -(_boosted_pim->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pim->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_pim->Vect().Unit();

  alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pim->Vect() < 0) alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

  //α[pπ+][p'π−]
  /// 2
  a_gamma = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_prot->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_prot->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_prot->Vect().Unit();

  alpha_PIpPIm_pipf = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_prot->Vect() < 0) alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf;
  //α[pp'][π+π−]

  /// 3
  a_gamma = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_pip->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pip->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()), 2)));
  b_beta = -(_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pim->Vect().Unit() + b_beta * _boosted_pip->Vect().Unit();

  alpha_PPIm_piPIp = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pip->Vect() < 0) alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp;

  _alpha_ppip_pipim = alpha_PPIp_piPIm;
  _alpha_pippim_pipf = alpha_PIpPIm_pipf;
  _alpha_ppim_pipip = alpha_PPIm_piPIp;
}

float Reaction::alpha_ppip_pipim() {  // pipim bhaneko proton initial  pim ho?
  if (_alpha_ppip_pipim != _alpha_ppip_pipim) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_ppip_pipim;
  else
    return NAN;
}
float Reaction::alpha_pippim_pipf() {  // alpha P (proton initial proton final)
  if (_alpha_pippim_pipf != _alpha_pippim_pipf) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_pippim_pipf;
  else
    return NAN;
}
float Reaction::alpha_ppim_pipip() {  // alpha pip (proton initial pip)
  if (_alpha_ppim_pipip != _alpha_ppim_pipip) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_ppim_pipip;
  else
    return NAN;
}
