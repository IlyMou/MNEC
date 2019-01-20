#ifndef _SCPACESCHEME_CPP
#include "SpaceScheme.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
SpaceScheme::SpaceScheme(DataFile* data_file) :
_norm_l2(data_file->Get_norm_l2()), _N(data_file->Get_N_mesh()),
_results(data_file->Get_results()), _data_file(data_file)
{
	system(("mkdir -p ./" + _results).c_str());
	_U.resize(5, _N);
	_F.resize(5, _N);
	_Ul.resize(5); Ul.resize(5);
	_Ur.resize(5); Ur.resize(5);
	check_order.setZero(5,2*_N);
	_stab =false;
}

// Construit la condition initiale au centre des triangles
MatrixXd SpaceScheme::InitialCondition()
{
	double dx = 1./_N;
	if(_data_file->Get_initial_condition_choice()=="riemann")
	{
		_x0 = _data_file->Get_param_x0();
		double hl = _data_file->Get_param_a();
		double ul = _data_file->Get_param_b();
		double vl = _data_file->Get_param_c();
		double al = _data_file->Get_param_d();
		double bl = _data_file->Get_param_e();
		double hr = _data_file->Get_param_f();
		double ur = _data_file->Get_param_g();
		double vr = _data_file->Get_param_h();
		double ar = _data_file->Get_param_i();
		double br = _data_file->Get_param_j();
		Ul << hl, ul, vl, al, bl;
		Ur << hr, ur, vr, ar, br;
		_Ul << hl, ul*hl, vl*hl, al*hl, bl*hl;
		_Ur << hr, ur*hr, vr*hr, ar*hr, br*hr;
		for (int i = 0 ; i < _N ; i++)
		{
			if((0.5+i)*dx<=_x0)
			{
				_U(0,i) = hl;
				_U(1,i) = ul*hl;
				_U(2,i) = vl*hl;
				_U(3,i) = al*hl;
				_U(4,i) = bl*hl;
			}
			else
			{
				_U(0,i) = hr;
				_U(1,i) = ur*hr;
				_U(2,i) = vr*hr;
				_U(3,i) = ar*hr;
				_U(4,i) = br*hr;
			}
		}
	}
	return _U;
}

double SpaceScheme::vp_b(const Eigen::MatrixXd& sol, int i)
{
	double b=0,bt=0;
	double ul, ur, al, ar, hl, hr;
	if ((i!=0)&&(i!=_N))
	{
		hl = sol(0,i-1); hr = sol(0,i);
		if((abs(sol(0,i))>1e-14)&&(abs(sol(0,i-1))>1e-14))
		{
			ul = sol(1,i-1)/sol(0,i-1); ur = sol(1,i)/sol(0,i);
			al = sol(3,i-1)/sol(0,i-1);	ar = sol(3,i)/sol(0,i);
		}
		if (abs(sol(0,i))<1e-14)
		{
			if(i<=_N/2) { ur = Ul(1);	ar = Ul(3); }
			else { ur = Ur(1);	ar = Ur(3); }
		}
		if (abs(sol(0,i-1))<1e-14)
		{
			if(i<=_N/2) { ul = Ul(1); al = Ul(3); }
			else { ul = Ur(1); al = Ur(3); } 
		}
	}
	else if (i==0)
	{
		hl = Ul(0); hr = sol(0,0);
		if(abs(hr)>1e-14)
		{
			ur = sol(1,0)/sol(0,0);
		  ar = sol(3,0)/sol(0,0);
		}
		else
		{
			ur = Ul(1);
			ar = Ul(3);
		}
		ul = Ul(1); al = Ul(3);
	}
	else
	{
		hl = sol(0,_N-1); hr = Ur(0);
		if(abs(hl)>1e-14)
		{
			ul = sol(1,_N-1)/sol(0,_N-1);
			al = sol(3,_N-1)/sol(0,_N-1);
		}
		else
		{
			ul = Ur(1);
			al = Ur(3);
		}
		ur = Ur(1); ar = Ur(3);
	}

	b  = max(abs(ul-sqrt(al*al+g*hl)),ul+sqrt(al*al+g*hl));
	bt = max(abs(ur-sqrt(ar*ar+g*hr)),ur+sqrt(ar*ar+g*hr));
	b  = max(b, bt);

	return b;
}


// -----------------------------------------
//   Rusanov ordre 1
// -----------------------------------------
Rusanov1::Rusanov1(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{}

// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
void Rusanov1::BuildF(const double& t, const Eigen::MatrixXd& sol)
{
	double bl = 0, br = 0; _bmax = 0;

	for(int i = 0; i<_N; i++)
	{
		_F.col(i) = 0.5*(Flux_R(sol,i)-Flux_L(sol,i));

		bl = vp_b(sol,i);
		br = vp_b(sol,i+1);
		if ((i!=0)&&(i!=_N-1))
		{
			_F.col(i) += 0.5*( bl*(sol.col(i)-sol.col(i-1)) - br*(sol.col(i+1)-sol.col(i)) );
		}
		else if (i==0)
		{
			_F.col(i) += 0.5*( bl*(sol.col(i)-_Ul) - br*(sol.col(i+1)-sol.col(i)) );
		}
		else
		{
			_F.col(i) += 0.5*( bl*(sol.col(i)-sol.col(i-1)) - br*(_Ur-sol.col(i)) );
		}

		if(_bmax<bl)
			_bmax = bl;
		if(_bmax<br)
			_bmax = br;
	}
	_F = -_N*_F;
}

VectorXd Rusanov1::Flux_R(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi;
	Fi.setZero(5);

	if (i!=_N-1)
	{
		if(abs(sol(0,i+1))>1e-14)
		{
			Fi(0) = sol(1,i+1);
			Fi(1) = ( sol(1,i+1)*sol(1,i+1) - sol(3,i+1)*sol(3,i+1) )/sol(0,i+1) + 0.5*g*sol(0,i+1)*sol(0,i+1) ;
			Fi(2) = ( sol(1,i+1)*sol(2,i+1) - sol(3,i+1)*sol(4,i+1) )/sol(0,i+1);
			Fi(4) = ( sol(1,i+1)*sol(4,i+1) - sol(2,i+1)*sol(3,i+1) )/sol(0,i+1);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = sol(3,i+1)*(sol(1,i)/sol(0,i));
			Fi(4) += sol(3,i+1)*(sol(2,i)/sol(0,i));
		}
	}
	else
	{
		Fi(0) = _Ur(1);
		Fi(1) = ( Ur(1)*Ur(1) - Ur(3)*Ur(3) )*Ur(0) + 0.5*g*_Ur(0)*_Ur(0);
		Fi(2) = ( Ur(1)*Ur(2) - Ur(3)*Ur(4) )*Ur(0);
		Fi(4) = ( Ur(1)*Ur(4) - Ur(2)*Ur(3) )*Ur(0);
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = _Ur(3)*(sol(1,i)/sol(0,i));
			Fi(4) += _Ur(3)*(sol(2,i)/sol(0,i));
		}
	}

	return Fi;
}

VectorXd Rusanov1::Flux_L(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi;
	Fi.setZero(5);

	if (i!=0)
	{
		if(abs(sol(0,i-1))>1e-14)
		{
			Fi(0) = sol(1,i-1);
			Fi(1) = ( sol(1,i-1)*sol(1,i-1) - sol(3,i-1)*sol(3,i-1) )/sol(0,i-1) + 0.5*g*sol(0,i-1)*sol(0,i-1);
			Fi(2) = ( sol(1,i-1)*sol(2,i-1) - sol(3,i-1)*sol(4,i-1) )/sol(0,i-1);
			Fi(4) = ( sol(1,i-1)*sol(4,i-1) - sol(2,i-1)*sol(3,i-1) )/sol(0,i-1);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = sol(3,i-1)*(sol(1,i)/sol(0,i));
			Fi(4) += sol(3,i-1)*(sol(2,i)/sol(0,i));
		}
	}
	else
	{
		Fi(0) = _Ul(1);
		Fi(1) = ( Ul(1)*Ul(1) - Ul(3)*Ul(3) )*Ul(0) + 0.5*g*_Ul(0)*_Ul(0);
		Fi(2) = ( Ul(1)*Ul(2) - Ul(3)*Ul(4) )*Ul(0);
		Fi(4) = ( Ul(1)*Ul(4) - Ul(2)*Ul(3) )*Ul(0);
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = _Ul(3)*(sol(1,i)/sol(0,i));
			Fi(4) += _Ul(3)*(sol(2,i)/sol(0,i));
		}
	}

	return Fi;
}


// -----------------------------------------
//   Rusanov ordre 2
// -----------------------------------------
Rusanov2::Rusanov2(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{
	if(data_file->Get_order() == 0)
	{
		_stab = true;
		_rO1 = new Rusanov1(data_file);
		_rO1->InitialCondition();
	}
}

// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
void Rusanov2::BuildF(const double& t, const Eigen::MatrixXd& sol)
{
	double bl = 0, br = 0; _bmax = 0;

	for(int i = 0; i<_N; i++)
	{
		_F.col(i) = 0.5*(Flux_R(sol,i)-Flux_L(sol,i));

		bl = vp_b(sol,i);
		br = vp_b(sol,i+1);
		if ((i!=0)&&(i<_N-2))
		{
			_F.col(i) += 0.5*( br*( sol.col(i+2)-2*sol.col(i+1)+sol.col(i) )
												-bl*( sol.col(i+1)-2*sol.col(i)+sol.col(i-1) ) );
		}
		else if (i==0)
		{
			_F.col(i) += 0.5*( br*( sol.col(i+2)-2*sol.col(i+1)+sol.col(i) )
												-bl*( sol.col(i+1)-2*sol.col(i)+_Ul ) );
		}
		else if (i==_N-2)
		{
			_F.col(i) += 0.5*( br*( _Ur-2*sol.col(i+1)+sol.col(i) )
												-bl*( sol.col(i+1)-2*sol.col(i)+sol.col(i-1) ) );
		}
		else
		{
			_F.col(i) += 0.5*( br*( -_Ur+sol.col(i) )
												-bl*( _Ur-2*sol.col(i)+sol.col(i-1) ) );
		}

		if(_bmax<bl)
			_bmax = bl;
		if(_bmax<br)
			_bmax = br;
	}
	_F = -_N*_F;


	if(_stab)
	{
		double phi;
		MatrixXd F;
		_rO1->BuildF(t, sol);
		F = _rO1->GetF();
		for(int j = 0; j<_N; j++)
		{
			for(int i = 0; i<5; i++){
				phi = limPente(sol,i,j);
				_F(i,j) = (1-phi)*F(i,j) + phi*_F(i,j);
				check_order(i,2*j) = phi;
				check_order(i,2*j+1) = phi;
			}
		}
	}
}

VectorXd Rusanov2::Flux_R(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi, solt;
	Fi.setZero(5);

	if (i<_N-2)
	{
		solt = 0.5*( 3*sol.col(i+1)-sol.col(i+2) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0) ;
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( sol.col(i+1)+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0) ;
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}
	else if (i==_N-2)
	{
		solt = 0.5*( 3*sol.col(i+1)-_Ur);
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( sol.col(i+1)+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}
	else
	{
		solt = _Ur;
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( _Ur+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);

		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}

	return Fi;
}

VectorXd Rusanov2::Flux_L(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi, solt;
	Fi.setZero(5);

	if ((i!=0)&&(i!=_N-1))
	{
		solt = 0.5*( 3*sol.col(i)-sol.col(i+1) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( sol.col(i-1)+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}
	else if (i==0)
	{
		solt = 0.5*( 3*sol.col(i)-sol.col(i+1) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3)  = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( _Ul+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}
	else
	{
		solt = 0.5*( 3*sol.col(i)-_Ur );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) = solt(1);
			Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}

		solt = 0.5*( sol.col(i-1)+sol.col(i) );
		if(abs(solt(0))>1e-14)
		{
			Fi(0) += solt(1);
			Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
			Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
			Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0);
		}
		if(abs(sol(0,i))>1e-14)
		{
			Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
			Fi(4) += solt(3)*(sol(2,i)/sol(0,i));
		}
	}

	return Fi;
}

double Rusanov2::limPente(const Eigen::MatrixXd& sol, int var, int i)
{
	double phi, theta;
	double u, v, w;

	v = sol(var,i);
	if((i > 0)&&(i<_N-1))
	{
		u = sol(var,i-1);
		w = sol(var,i+1);
	}
	else if (i == 0)
	{
		u = _Ul(var);
		w = sol(var,i+1);
	}
	else
	{
		u = sol(var,i-1);
		w = _Ur(var);
	}

	if((abs(v-u)<1e-12)&&(abs(w-v)<1e-12))
		theta = 1;
	else if(abs(w-v)<2.5e-13)
		theta = 0;
	else
		theta = (v-u)/(w-u);

	if((theta<0.25)||(theta >4))
		theta = 0;

	phi = (theta + abs(theta))/(1.+theta);

	return phi;
}


// -----------------------------------------
//   Relaxation method 1
// -----------------------------------------
WRS1::WRS1(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{
	_dt = data_file->Get_dt();
	_ISl.setZero(9,6);
	_ISr.setZero(9,6);
	_stab = false; _o2 = false;
	if(data_file->Get_order() == 2)
		_o2 = true;
	if(data_file->Get_order() == 0){
		_o2 = true; _stab = true;
	}
}

VectorXd WRS1::vp_c(const Eigen::MatrixXd& IS)
{
	double hr_, hl_, hr, hl;
	double cl, cr, cal, car;
	double ur, ul, u_;

	hl = IS(0,0); hr = IS(0,5);
	hl_= IS(0,1); hr_= IS(0,4);
	ul = IS(1,0); ur = IS(1,5);	u_ = IS(1,3);
	cl = IS(7,0); cr = IS(7,5);
	cal= IS(8,0); car= IS(8,5);

	VectorXd vp(5);
	vp(0) = ul; vp(1) = u_;
	vp(2) = u_;
	vp(3) = u_;	vp(4) = ur;

	if(hl>1e-12)
		vp(0) -= cl/hl;
	if(hl_>1e-12)
		vp(1) -= cal/hl_;
	if(hr_>1e-12)
		vp(3) += car/hr_;
	if(hr>1e-12)
		vp(4) += cr/hr;

	return vp;
}

void WRS1::BuildF(const double& t, const Eigen::MatrixXd& sol)
{
	double bl = 0, br = 0; _bmax = 0;
	_Sigl.setZero(5);
 	_Sigr.setZero(5);

	_ISl = interState(sol, 0);
	for(int i = 0; i<_N; i++)
	{
		if(i!=0) _ISl = _ISr;
		_ISr = interState(sol, i+1);

		_Sigl = vp_c(_ISl);	bl = max(abs(_Sigl(4)),abs(_Sigl(0)));
		_Sigr = vp_c(_ISr);	br = max(abs(_Sigr(4)),abs(_Sigr(0)));

		if(_bmax<bl)
			_bmax = bl;
		if(_bmax<br)
			_bmax = br;

		_F.col(i) = (Flux_R(sol,i)-Flux_L(sol,i));
	}

	//cout << "bmax = " << _bmax << endl;

	_F = -_N*_F;
}

VectorXd WRS1::Flux_R(const Eigen::MatrixXd& sol, int i)
{
	double dx = 1./_N;
	VectorXd Fi, vp = _Sigr;
	Fi.setZero(5);

	if (i!=_N-1)
	{
		if(abs(sol(0,i))>1e-14)
		{
			Fi(0) = sol(1,i);
			Fi(1) = ( sol(1,i)*sol(1,i) - sol(3,i)*sol(3,i) )/sol(0,i) + 0.5*g*sol(0,i)*sol(0,i) ;
			Fi(2) = ( sol(1,i)*sol(2,i) - sol(3,i)*sol(4,i) )/sol(0,i);
			Fi(4) = ( sol(1,i)*sol(4,i) - sol(2,i)*sol(3,i) )/sol(0,i);
		}
	}
	else
	{
		Fi(0) = _Ur(1);
		Fi(1) = ( Ur(1)*Ur(1) - Ur(3)*Ur(3) )*Ur(0) + 0.5*g*_Ur(0)*_Ur(0);
		Fi(2) = ( Ur(1)*Ur(2) - Ur(3)*Ur(4) )*Ur(0);
		Fi(4) = ( Ur(1)*Ur(4) - Ur(2)*Ur(3) )*Ur(0);
	}



	if( _N*_dt*vp(4) < -0.5 )
  {
    cout << endl << "The CFL is not respected ! " << endl; abort();
  }

	MatrixXd IS = _ISr;
	for(int k = 1; k<5; k++)
	{
		for(int j = 0; j<6; j++)
		{
			IS(k,j) = IS(k,j)*IS(0,j);
		}
	}

	if( vp(4)<0 )
	{
		Fi(0) -= -vp(4)*IS(0,5) + (vp(4)-vp(2))*IS(0,4) + (vp(2)-vp(0))*IS(0,2)
						+ vp(0)*IS(0,0);
		Fi(1) -= -vp(4)*IS(1,5) + (vp(4)-vp(2))*IS(1,4) + (vp(2)-vp(0))*IS(1,1)
						+ vp(0)*IS(1,0);
		Fi(2) -= -vp(4)*IS(2,5) + (vp(4)-vp(3))*IS(2,4) + (vp(3)-vp(2))*IS(2,3)
						+ (vp(2)-vp(1))*IS(2,2) + (vp(1)-vp(0))*IS(2,1) + vp(0)*IS(2,0);
		Fi(3) -= -vp(4)*IS(3,5) + (vp(4)-vp(2))*IS(3,4) + (vp(2)-vp(0))*IS(3,2)
						+ vp(0)*IS(3,0);
		Fi(4) -= -vp(4)*IS(4,5) + (vp(4)-vp(3))*IS(4,4) + (vp(3)-vp(2))*IS(4,3)
						+ (vp(2)-vp(1))*IS(4,2) + (vp(1)-vp(0))*IS(4,1) + vp(0)*IS(4,0);
	}
	else if( vp(3) < 0)
	{
		Fi(0) -= -vp(2)*IS(0,4) + (vp(2)-vp(0))*IS(0,2) + vp(0)*IS(0,0);
		Fi(1) -= -vp(2)*IS(1,4) + (vp(2)-vp(0))*IS(1,1) + vp(0)*IS(1,0);
		Fi(2) -= -vp(3)*IS(2,4) + (vp(3)-vp(2))*IS(2,3) + (vp(2)-vp(1))*IS(2,2)
						+ (vp(1)-vp(0))*IS(2,1) + vp(0)*IS(2,0);
		Fi(3) -= -vp(2)*IS(3,4) + (vp(2)-vp(0))*IS(3,2) + vp(0)*IS(3,0);
		Fi(4) -= -vp(3)*IS(4,4) + (vp(3)-vp(2))*IS(4,3) + (vp(2)-vp(1))*IS(4,2)
						+ (vp(1)-vp(0))*IS(4,1) + vp(0)*IS(4,0);
	}
	else if( vp(2) < 0)
	{
		Fi(0) -= -vp(2)*IS(0,4) + (vp(2)-vp(0))*IS(0,2) + vp(0)*IS(0,0);
		Fi(1) -= -vp(2)*IS(1,4) + (vp(2)-vp(0))*IS(1,1) + vp(0)*IS(1,0);
		Fi(2) -= -vp(2)*IS(2,3) + (vp(2)-vp(1))*IS(2,2) + (vp(1)-vp(0))*IS(2,1)
						+ vp(0)*IS(2,0);
		Fi(3) -= -vp(2)*IS(3,4) + (vp(2)-vp(0))*IS(3,2) + vp(0)*IS(3,0);
		Fi(4) -= -vp(2)*IS(4,3) + (vp(2)-vp(1))*IS(4,2) + (vp(1)-vp(0))*IS(4,1)
						+ vp(0)*IS(4,0);
	}
	else if( vp(1) < 0)
	{
		Fi(0) -= -vp(0)*IS(0,2) + vp(0)*IS(0,0);
		Fi(1) -= -vp(0)*IS(1,1) + vp(0)*IS(1,0);
		Fi(2) -= -vp(1)*IS(2,2) + (vp(1)-vp(0))*IS(2,1) + vp(0)*IS(2,0);
		Fi(3) -= -vp(0)*IS(3,2) + vp(0)*IS(3,0);
		Fi(4) -= -vp(1)*IS(4,2) + (vp(1)-vp(0))*IS(4,1) + vp(0)*IS(4,0);
	}
	else if( vp(0) < 0)
	{
		Fi(0) -= -vp(0)*IS(0,2) + vp(0)*IS(0,0);
		Fi(1) -= -vp(0)*IS(1,1) + vp(0)*IS(1,0);
		Fi(2) -= -vp(0)*IS(2,1) + vp(0)*IS(2,0);
		Fi(3) -= -vp(0)*IS(3,2) + vp(0)*IS(3,0);
		Fi(4) -= -vp(0)*IS(4,1) + vp(0)*IS(4,0);
	}

	return Fi;
}

VectorXd WRS1::Flux_L(const Eigen::MatrixXd& sol, int i)
{
	double dx = 1./_N;
	VectorXd Fi, vp = _Sigl;
	Fi.setZero(5);

	if (i!=0)
	{
		if(abs(sol(0,i))>1e-14)
		{
			Fi(0) = sol(1,i);
			Fi(1) = ( sol(1,i)*sol(1,i) - sol(3,i)*sol(3,i) )/sol(0,i) + 0.5*g*sol(0,i)*sol(0,i);
			Fi(2) = ( sol(1,i)*sol(2,i) - sol(3,i)*sol(4,i) )/sol(0,i);
			Fi(4) = ( sol(1,i)*sol(4,i) - sol(2,i)*sol(3,i) )/sol(0,i);
		}
	}
	else
	{
		Fi(0) = _Ul(1);
		Fi(1) = ( Ul(1)*Ul(1) - Ul(3)*Ul(3) )*Ul(0) + 0.5*g*_Ul(0)*_Ul(0);
		Fi(2) = ( Ul(1)*Ul(2) - Ul(3)*Ul(4) )*Ul(0);
		Fi(4) = ( Ul(1)*Ul(4) - Ul(2)*Ul(3) )*Ul(0);
	}


	if( _N*_dt*vp(4) > 0.5 )
  {
    cout << endl << "The CFL is not respected ! " << endl; abort();
  }

	MatrixXd IS = _ISl;
	for(int k = 1; k<5; k++)
	{
		for(int j = 0; j<6; j++)
		{
			IS(k,j) = IS(k,j)*IS(0,j);
		}
	}


	if( vp(0) > 0 )
	{
		Fi(0) += vp(0)*IS(0,0) + (vp(2)-vp(0))*IS(0,1) + (vp(4)-vp(2))*IS(0,3)
						-vp(4)*IS(0,5);
		Fi(1) += vp(0)*IS(1,0) + (vp(2)-vp(0))*IS(1,2) + (vp(4)-vp(2))*IS(1,3)
		 				-vp(4)*IS(1,5);
		Fi(2) += vp(0)*IS(2,0) + (vp(1)-vp(0))*IS(2,1) + (vp(2)-vp(1))*IS(2,2)
		 				+(vp(3)-vp(2))*IS(2,3) + (vp(4)-vp(3))*IS(2,4) - vp(4)*IS(2,5);
		Fi(3) += vp(0)*IS(3,0) + (vp(2)-vp(0))*IS(3,1) + (vp(4)-vp(2))*IS(3,3)
						-vp(4)*IS(3,5);
		Fi(4) += vp(0)*IS(4,0) + (vp(1)-vp(0))*IS(4,1) + (vp(2)-vp(1))*IS(4,2)
		 				+(vp(3)-vp(2))*IS(4,3) + (vp(4)-vp(3))*IS(4,4) - vp(4)*IS(4,5);
	}
	else if( vp(1) > 0)
	{
		Fi(0) += vp(2)*IS(0,1) + (vp(4)-vp(2))*IS(0,3) - vp(4)*IS(0,5);
		Fi(1) += vp(2)*IS(1,2) + (vp(4)-vp(2))*IS(1,3) - vp(4)*IS(1,5);
		Fi(2) += vp(1)*IS(2,1) + (vp(2)-vp(1))*IS(2,2) + (vp(3)-vp(2))*IS(2,3)
						+(vp(4)-vp(3))*IS(2,4) - vp(4)*IS(2,5);
		Fi(3) += vp(2)*IS(3,1) + (vp(4)-vp(2))*IS(3,3) - vp(4)*IS(3,5);
		Fi(4) += vp(1)*IS(4,1) + (vp(2)-vp(1))*IS(4,2) + (vp(3)-vp(2))*IS(4,3)
						+(vp(4)-vp(3))*IS(4,4) - vp(4)*IS(4,5);
	}
	else if( vp(2) > 0)
	{
		Fi(0) += vp(2)*IS(0,1) + (vp(4)-vp(2))*IS(0,3) - vp(4)*IS(0,5);
		Fi(1) += vp(2)*IS(1,3) + (vp(4)-vp(2))*IS(1,3) - vp(4)*IS(1,5);
		Fi(2) += vp(2)*IS(2,2) + (vp(3)-vp(2))*IS(2,3) + (vp(4)-vp(3))*IS(2,4)
						-vp(4)*IS(2,5);
		Fi(3) += vp(2)*IS(3,1) + (vp(4)-vp(2))*IS(3,3) - vp(4)*IS(3,5);
		Fi(4) += vp(2)*IS(4,2) + (vp(3)-vp(2))*IS(4,3) + (vp(4)-vp(3))*IS(4,4)
						-vp(4)*IS(4,5);
	}
	else if( vp(3) > 0)
	{
		Fi(0) += vp(4)*IS(0,3) - vp(4)*IS(0,5);
		Fi(1) += vp(4)*IS(1,3) - vp(4)*IS(1,5);
		Fi(2) += vp(3)*IS(2,3) + (vp(4)-vp(3))*IS(2,4) - vp(4)*IS(2,5);
		Fi(3) += vp(4)*IS(3,3) - vp(4)*IS(3,5);
		Fi(4) += vp(3)*IS(4,3) + (vp(4)-vp(3))*IS(4,4) - vp(4)*IS(4,5);
	}
	else if( vp(4) > 0)
	{
		Fi(0) += vp(4)*IS(0,3) - vp(4)*IS(0,5);
		Fi(1) += vp(4)*IS(1,3) - vp(4)*IS(1,5);
		Fi(2) += vp(4)*IS(2,4) - vp(4)*IS(2,5);
		Fi(3) += vp(4)*IS(3,3) - vp(4)*IS(3,5);
		Fi(4) += vp(4)*IS(4,4) - vp(4)*IS(4,5);
	}

	return Fi;
}

MatrixXd WRS1::interState(const Eigen::MatrixXd& sol, int i)
{
	double hl, hr, hl_, hr_;
	double ul, ur, u_;
	double vl, vr, v_;
	double al, ar, al_, ar_;
	double bl, br, bl_, br_;
	double pl, pr, p_;
	double ptl, ptr, pt_;
	double cl, cr, cl_, cr_;
	double cal, car, cal_, car_;
	double sl, sr;

	MatrixXd IS;
	IS.setZero(9,6);

	//	Definition des variables gauche et droite
	if((i>0)&&(i<_N))
	{
		hl = sol(0,i-1);
		hr = sol(0,i);
		if((abs(hl)>1e-12)&&(abs(hr)>1e-12))
		{
			ul = sol(1,i-1)/hl; ur = sol(1,i)/hr;
			vl = sol(2,i-1)/hl; vr = sol(2,i)/hr;
			al = sol(3,i-1)/hl;	ar = sol(3,i)/hr;
			bl = sol(4,i-1)/hl;	br = sol(4,i)/hr;
		}
		else if ((abs(hl)<1e-12)&&(abs(hr)>1e-12))
		{
			ul = Ul(1); ur = sol(1,i)/hr;
			vl = Ul(2); vr = sol(2,i)/hr;
			al = Ul(3); ar = sol(3,i)/hr;
			bl = Ul(4); br = sol(4,i)/hr;
		}
		else if ((abs(hl)>1e-12)&&(abs(hr)<1e-12))
		{
			ul = sol(1,i-1)/hl; ur = Ur(1);
			vl = sol(2,i-1)/hl; vr = Ur(2);
			al = sol(3,i-1)/hl;	ar = Ur(3);
			bl = sol(4,i-1)/hl;	br = Ur(4);
		}
		else
		{
			if(i>=_N/2)
			{
				ul = Ur(1); ur = Ur(1);
				ul = Ur(2); ur = Ur(2);
				al = Ur(3);	ar = Ur(3);
				bl = Ur(4);	br = Ur(4);
			}
			else
			{
				ul = Ul(1); ur = Ul(1);
				ul = Ul(2); ur = Ul(2);
				al = Ul(3);	ar = Ul(3);
				bl = Ul(4);	br = Ul(4);
			}
		}
	}
	else if(i==0)
	{
		hr = sol(0,i); hl = _Ul(0);
		ul = Ul(1); vl = Ul(2);
		al = Ul(3); bl = Ul(4);
		if(hr>1e-12)
		{
			ur = sol(1,i)/hr; vr = sol(2,i)/hr;
			ar = sol(3,i)/hr; br = sol(4,i)/hr;
		}
		else
		{
			ur = Ur(1); vr = Ur(2);
			ar = Ur(3);	br = Ur(4);
		}
	}
	else
	{
		hl = sol(0,i-1); hr = Ur(0);
		ur = Ur(1); vr = Ur(2);
		ar = Ur(3); br = Ur(4);
		if(hl>1e-12)
		{
			ul = sol(1,i-1)/hl; vl = sol(2,i-1)/hl;
			al = sol(3,i-1)/hl; bl = sol(4,i-1)/hl;
		}
		else
		{
			ul = Ul(1); vl = Ul(2);
			al = Ul(3); bl = Ul(4);
		}
	}

	sl = sqrt(al*al+g*hl);
	sr = sqrt(ar*ar+g*hr);
	pl = (g*hl*0.5 - al*al)*hl;
	pr = (g*hr*0.5 - ar*ar)*hr;
	ptl = -hl*al*bl;
	ptr = -hr*ar*br;

	cl = 0; cr = 0;
	if(hl>1e-12)
		cl = hl*(sl + 1.5*(max(0.,ul-ur) + max(0.,pr - pl)/(hl*sl + hr*sr)));
	if(hr>1e-12)
		cr = hr*(sr + 1.5*(max(0.,ul-ur) + max(0.,pl - pr)/(hl*sl + hr*sr)));

	cal = hl*abs(al);
	car = hr*abs(ar);


	// Calcul des etats intermediaires
	hl_ = 0; hr_ = 0;
	if(hl>1e-12)
		hl_ = hl/( 1. + hl*( cr*(ur-ul) + pl-pr )/( cl*(cl+cr) ) );
	if(hr>1e-12)
		hr_ = hr/( 1. + hr*( cl*(ur-ul) + pr-pl )/( cr*(cl+cr) ) );

	u_  = 0; v_  = 0;
	p_  = 0; pt_ = 0;
	bl_ = bl; br_ = br;
	if((hr>1e-12)||(hl>1e-12))
	{
		u_ = ( cl*ul + cr*ur + pl - pr )/( cl + cr );
		v_ = ( cal*vl + car*vr + ptl - ptr )/( cal + car );

		p_ = ( cr*pl + cl*pr - cl*cr*(ur-ul) )/( cl+cr );
		pt_= ( car*ptl + cal*ptr - cal*car*(vr-vl) )/( cal+car );

		if(abs(car)>1e-12)
			br_= br + (ptr - ptl + cal*(vr-vl)) * hr*ar/(car*(cal+car));
		if(abs(cal)>1e-12)
			bl_= bl + (ptl - ptr + car*(vr-vl)) * hl*al/(cal*(cal+car));
	}

	al_ = al; ar_ = ar;
	if(hl>1e-12)
		al_ = al*hl/hl_;
 	if(hr>1e-12)
		ar_ = ar*hr/hr_;

	// Faux
	cl_ = cl; cr_ = cr;
	cal_ = cal; car_ = car;

	// Mise en mémoire des variables
	IS << hl, hl_ , hl_ , hr_ , hr_ , hr ,
				ul , u_  , u_  , u_  , u_  , ur ,
				vl , vl  , v_  , v_  , vr  , vr ,
				al , al_ , al_ , ar_ , ar_ , ar ,
				bl , bl  , bl_ , br_ , br  , br ,
				pl , p_  , p_  , p_  , p_  , pr ,
				ptl, ptl , pt_ , pt_ , ptr , ptr,
				cl , cl_ , cl_ , cr_ , cr_ , cr ,
				cal, cal_, cal_, car_, car_, car;

	// if((i == _N/2-1)||(i == _N/2)||(i == _N/2+1))
	// 	cout << endl << endl << " i = " << i << " " << endl << IS << endl;
	return IS;
}


// -----------------------------------------
//   Relaxation method
// -----------------------------------------
WRS2::WRS2(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{

}

// Une étape du schéma en temps
void WRS2::BuildF(const double& t, const Eigen::MatrixXd& sol)
{

}

Eigen::VectorXd WRS2::Flux_R(const Eigen::MatrixXd& sol, int i)
{

}

Eigen::VectorXd WRS2::Flux_L(const Eigen::MatrixXd& sol, int i)
{

}

Eigen::VectorXd WRS2::vp_c(const Eigen::MatrixXd& IS)
{

}

Eigen::MatrixXd WRS2::interState(const Eigen::MatrixXd& sol, int i)
{

}


// -----------------------------------------
//   Solution
// -----------------------------------------
// Sauvegarde la solution
void SpaceScheme::SaveSol(const Eigen::MatrixXd& sol, int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".dat";

	assert((sol.cols() == _N) && "The size of the solution matrix is not the same than the number of meshes !");
	assert((sol.rows() == 5) && "The size of the solution matrix is not the same than the number of variables !");

	double dx = 1./_N;

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(12);

	for (size_t i = 0; i < _N; i++)
	{
		if(abs(sol(0,i))>1e-14)
		{
			solution << i*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																					<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
			solution << (i+1)*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																							<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
		}
		else
		{
			if((i+0.5)*dx>_x0)
			{
				solution << i*dx << " " << sol(0,i) << " " << Ur(1) << " " << Ur(2)
																						<< " " << Ur(3) << " " << Ur(4) << endl;
				solution << (i+1)*dx << " " << sol(0,i) << " " << Ur(1) << " " << Ur(2)
																								<< " " << Ur(3) << " " << Ur(4) << endl;
			}
			else
			{
				solution << i*dx << " " << sol(0,i) << " " << Ul(1) << " " << Ul(2)
																						<< " " << Ul(3) << " " << Ul(4) << endl;
				solution << (i+1)*dx << " " << sol(0,i) << " " << Ul(1) << " " << Ul(2)
																								<< " " << Ul(3) << " " << Ul(4) << endl;
			}
		}
	}

  solution << endl;
	solution.close();
}

// Sauvegarde la derniere solution
void SpaceScheme::SaveSol(const Eigen::MatrixXd& sol)
{
	string name_file = _results + "/sol.dat";

	assert((sol.cols() == _N) && "The size of the solution matrix is not the same than the number of meshes !");
	assert((sol.rows() == 5) && "The size of the solution matrix is not the same than the number of variables !");

	double dx = 1./_N;

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(12);

	solution << 0. << " " << Ul(0) << " " << Ul(1) << " " << Ul(2)
																	<< " " << Ul(3) << " " << Ul(4) << endl;
	for (size_t i = 0; i < _N; i++)
		if(abs(sol(0,i))>1e-14)
			solution << (i+0.5)*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																								<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
		else
			if((i+0.5)*dx>_x0)
				solution << (i+0.5)*dx << " " << sol(0,i) << " " << Ur(1) << " " << Ur(2)
																									<< " " << Ur(3) << " " << Ur(4) << endl;
			else
				solution << (i+0.5)*dx << " " << sol(0,i) << " " << Ul(1) << " " << Ul(2)
																									<< " " << Ul(3) << " " << Ul(4) << endl;

  solution << 1. << " " << Ur(0) << " " << Ur(1) << " " << Ur(2)
																	<< " " << Ur(3) << " " << Ur(4);

	solution.close();

	if(_stab)
	{
		name_file = _results + "/limPente.dat";
		ofstream order;
		order.open(name_file, ios::out);
		order.precision(12);
		for (int i = 0; i < _N; i++)
		{
			order << i*dx;
			for (int j = 0; j < 5; j++)
				order << " " << check_order(j,2*i);
			order << endl;

			order << (i+1)*dx;
			for (int j = 0; j < 5; j++)
				order << " " << check_order(j,2*i+1);
			order << endl;
		}
		order.close();
	}
}

// Solution exacte au centre des triangles
MatrixXd SpaceScheme::ExactSolution()
{
	int N = _data_file->Get_N_solE();
	double x,h,u,v,a,b;
	string name_file = _data_file->Get_file_solEx();

	ifstream solutionE;
	solutionE.open(name_file.data());
	if (!solutionE.is_open())
	{
		cout << "Unable to open file " << name_file << endl;
		abort();
	}

	MatrixXd exact_sol(5,N);
	solutionE >> x >> h >> u >> v >> a >> b;
	for (int i = 0 ; i < N ; i++)
	{
		solutionE >> x >> h >> u >> v >> a >> b;
		exact_sol(0,i) = h;
		exact_sol(1,i) = u;
		exact_sol(2,i) = v;
		exact_sol(3,i) = a;
		exact_sol(4,i) = b;
	}

	solutionE.close();

	return exact_sol;
}

void SpaceScheme::ComputeError(const Eigen::MatrixXd& sol)
{
	int N = _data_file->Get_N_solE(), k;
	double dx = 1./N, _dx = 1./_N;
	VectorXd error(5), eL2(5);
	MatrixXd exact_sol;
	exact_sol = ExactSolution();

	// --- Verification ---
	string name_file = _results + "/bruh.dat";
	ofstream bruh;
	bruh.open(name_file, ios::out);
	bruh.precision(12);


	bruh << 0. << " " << Ul(0) << " " << Ul(1) << " " << Ul(2)
															<< " " << Ul(3) << " " << Ul(4) << endl;
	for (size_t i = 0; i < N; i++)
		bruh << (i+0.5)*dx << " " << exact_sol(0,i) << " " << exact_sol(1,i) << " " << exact_sol(2,i)
																								<< " " << exact_sol(3,i) << " " << exact_sol(4,i) << endl;
  bruh << 1. << " " << Ur(0) << " " << Ur(1) << " " << Ur(2)
															<< " " << Ur(3) << " " << Ur(4);
	bruh.close();
	// --------------------

	error.setZero(5);; eL2.setZero(5);
	k=0;
	for (int i = 0 ; i < N ; i++)
	{
		if((i+1)*dx <= (k+1)*_dx)
		{
			error(0) += dx*(exact_sol(0,i) - sol(0,k))*(exact_sol(0,i) - sol(0,k));
			if(abs(sol(0,k))>1e-14)
				for(int j = 1; j<5; j++)
					error(j) += dx*(exact_sol(j,i) - sol(j,k)/sol(0,k))*(exact_sol(j,i) - sol(j,k)/sol(0,k));
			else
				for(int j = 1; j<5; j++)
					error(j) += dx*(exact_sol(j,i))*(exact_sol(j,i));
		}
		else
		{
			error(0) += ((k+1)*_dx  - i*dx)*(exact_sol(0,i) - sol(0,k))*(exact_sol(0,i) - sol(0,k));
			error(0) += ((i+1)*dx - (k+1)*_dx)*(exact_sol(0,i) - sol(0,k+1))*(exact_sol(0,i) - sol(0,k+1));

			if(abs(sol(0,k))>1e-14)
				for(int j = 1; j<5; j++)
					error(j) += ((k+1)*_dx  - i*dx)*(exact_sol(j,i) - sol(j,k)/sol(0,k))*(exact_sol(j,i) - sol(j,k)/sol(0,k));
			else
				for(int j = 1; j<5; j++)
					error(j) += ((k+1)*_dx  - i*dx)*(exact_sol(j,i))*(exact_sol(j,i));

			if(abs(sol(0,k+1))>1e-14)
				for(int j = 1; j<5; j++)
					error(j) += ((i+1)*dx - (k+1)*_dx)*(exact_sol(j,i) - sol(j,k+1)/sol(0,k+1))*(exact_sol(j,i) - sol(j,k+1)/sol(0,k+1));
			else
				for(int j = 1; j<5; j++)
					error(j) += ((i+1)*dx - (k+1)*_dx)*(exact_sol(j,i))*(exact_sol(j,i));
			k++;
		}

		for(int j=0; j<5; j++)
			eL2(j) += dx*exact_sol(j,i)*exact_sol(j,i);
	}

	cout << "Error L2 : " << sqrt(error(0) + error(1) + error(2) + error(3) + error(4))
													/sqrt(eL2(0) + eL2(1) + eL2(2) + eL2(3) + eL2(4)) << endl;
	cout << "  -- h = " << sqrt(error(0))/sqrt(eL2(0)) << endl;
	cout << "  -- u = " << sqrt(error(1))/sqrt(eL2(1)) << endl;
	cout << "  -- v = " << sqrt(error(2))/sqrt(eL2(2)) << endl;
	cout << "  -- a = " << sqrt(error(3))/sqrt(eL2(3)) << endl;
	cout << "  -- b = " << sqrt(error(4))/sqrt(eL2(4)) << endl;
}

#define _SCPACESCHEME_CPP
#endif
