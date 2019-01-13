#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme() : _pb(0)
{}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(DataFile* data_file, SpaceScheme* pb)
{
  _dt = data_file->Get_dt();
  _N = data_file->Get_N_mesh();
  _t = data_file->Get_t0();
  _pb = pb;
  _sol0 = pb->InitialCondition();
  _sol = _sol0;
}

void TimeScheme::SaveSolution(int n)
{
  _pb->SaveSol(_sol, n);
}

void TimeScheme::SaveSolution()
{
  _pb->SaveSol(_sol);
}

// Renvoie _sol (pratique pour vérifier la résolution)
const MatrixXd & TimeScheme::GetIterateSolution() const
{
  return _sol;
}

// Euler Explicite
void EulerScheme::Advance()
{
  _pb->BuildF(_t, _sol);
  _sol += _dt*_pb->GetF();
  _t += _dt;
}

// RungeKutta
void RK2::Advance()
{
  MatrixXd k1, k2;
  _pb->BuildF(_t, _sol);
  k1 = _sol + 0.5*_dt*_pb->GetF();
  _pb->BuildF(_t+0.5*_dt, k1);
  k2 = _pb->GetF();
  _sol += _dt*k2;
  _t += _dt;
}

// RungeKutta
void SSPRK2::Advance()
{
  MatrixXd k1, k2;
  _pb->BuildF(_t, _sol);
  k1 = _sol + _dt*_pb->GetF();
  _pb->BuildF(_t, k1);
  k2 = _pb->GetF();
  _sol = 0.5*_sol +0.5*k1 + 0.5*_dt*k2;
  _t += _dt;
}

void RK4::Advance()
{
  MatrixXd k1, k2, k3, k4;
  _pb->BuildF(_t, _sol);
  k1 = _pb->GetF();
  _pb->BuildF(_t+_dt/2., _sol+_dt/2.*k1);
  k2 = _pb->GetF();
  _pb->BuildF(_t+_dt/2., _sol+_dt/2.*k2);
  k3 = _pb->GetF();
  _pb->BuildF(_t+_dt, _sol+_dt*k3);
  k4 = _pb->GetF();
  _sol += _dt/6.*(k1 + 2.*k2 + 2.*k3 + k4);
  _t += _dt;
}


#define _TIME_SCHEME_CPP
#endif
