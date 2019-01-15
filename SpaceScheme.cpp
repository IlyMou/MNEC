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
	_Ul.resize(5);
	_Ur.resize(5);
}

// Construit la condition initiale au centre des triangles
MatrixXd SpaceScheme::InitialCondition()
{
	double dx = 1./_N;
	if(_data_file->Get_initial_condition_choice()=="riemann")
	{
		double x0 = _data_file->Get_param_x0();
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
		_Ul << hl, ul*hl, vl*hl, al*hl, bl*hl;
		_Ur << hr, ur*hr, vr*hr, ar*hr, br*hr;
		for (int i = 0 ; i < _N ; i++)
		{
			if((0.5+i)*dx<=x0)
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
		ul = sol(1,i-1)/sol(0,i-1); ur = sol(1,i)/sol(0,i);
		al = sol(3,i-1)/sol(0,i-1);	ar = sol(3,i)/sol(0,i);
		hl = sol(0,i-1);            hr = sol(0,i);
	}
	else if (i==0)
	{
		ul = _Ul(1)/_Ul(0); ur = sol(1,0)/sol(0,0);
		al = _Ul(3)/_Ul(0);	ar = sol(3,0)/sol(0,0);
		hl = _Ul(0);        hr = sol(0,0);
	}
	else
	{
		ul = sol(1,_N-1)/sol(0,_N-1); ur = _Ur(1)/_Ur(0);
		al = sol(3,_N-1)/sol(0,_N-1);	ar = _Ur(3)/_Ur(0);
		hl = sol(0,_N-1);             hr = _Ur(0);
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
	Fi.resize(5);

	if (i!=_N-1)
	{
		Fi(0) = sol(1,i+1);
		Fi(1) = ( sol(1,i+1)*sol(1,i+1) - sol(3,i+1)*sol(3,i+1) )/sol(0,i+1) + 0.5*g*sol(0,i+1)*sol(0,i+1) ;
		Fi(2) = ( sol(1,i+1)*sol(2,i+1) - sol(3,i+1)*sol(4,i+1) )/sol(0,i+1);
		Fi(3) = sol(3,i+1)*(sol(1,i)/sol(0,i));
		Fi(4) = ( sol(1,i+1)*sol(4,i+1) - sol(2,i+1)*sol(3,i+1) )/sol(0,i+1) + sol(3,i+1)*(sol(2,i)/sol(0,i));
	}
	else
	{
		Fi(0) = _Ur(1);
		Fi(1) = ( _Ur(1)*_Ur(1) - _Ur(3)*_Ur(3) )/_Ur(0) + 0.5*g*_Ur(0)*_Ur(0);
		Fi(2) = ( _Ur(1)*_Ur(2) - _Ur(3)*_Ur(4) )/_Ur(0);
		Fi(3) = _Ur(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( _Ur(1)*_Ur(4) - _Ur(2)*_Ur(3) )/_Ur(0) + _Ur(3)*(sol(2,i)/sol(0,i));
	}

	return Fi;
}

VectorXd Rusanov1::Flux_L(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi;
	Fi.resize(5);

	if (i!=0)
	{
		Fi(0) = sol(1,i-1);
		Fi(1) = ( sol(1,i-1)*sol(1,i-1) - sol(3,i-1)*sol(3,i-1) )/sol(0,i-1) + 0.5*g*sol(0,i-1)*sol(0,i-1);
		Fi(2) = ( sol(1,i-1)*sol(2,i-1) - sol(3,i-1)*sol(4,i-1) )/sol(0,i-1);
		Fi(3) = sol(3,i-1)*(sol(1,i)/sol(0,i));
		Fi(4) = ( sol(1,i-1)*sol(4,i-1) - sol(2,i-1)*sol(3,i-1) )/sol(0,i-1) + sol(3,i-1)*(sol(2,i)/sol(0,i));
	}
	else
	{
		Fi(0) = _Ul(1);
		Fi(1) = ( _Ul(1)*_Ul(1) - _Ul(3)*_Ul(3) )/_Ul(0) + 0.5*g*_Ul(0)*_Ul(0);
		Fi(2) = ( _Ul(1)*_Ul(2) - _Ul(3)*_Ul(4) )/_Ul(0);
		Fi(3) = _Ul(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( _Ul(1)*_Ul(4) - _Ul(2)*_Ul(3) )/_Ul(0) + _Ul(3)*(sol(2,i)/sol(0,i));
	}

	return Fi;
}


// -----------------------------------------
//   Rusanov ordre 1
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
			}
		}
	}
}

VectorXd Rusanov2::Flux_R(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi, solt;
	Fi.resize(5);

	if (i<_N-2)
	{
		solt = 0.5*( 3*sol.col(i+1)-sol.col(i+2) );

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0) ;
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( sol.col(i+1)+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0) ;
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
	}
	else if (i==_N-2)
	{
		solt = 0.5*( 3*sol.col(i+1)-_Ur);

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( sol.col(i+1)+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
	}
	else
	{
		solt = _Ur;

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( _Ur+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
	}

	return Fi;
}

VectorXd Rusanov2::Flux_L(const Eigen::MatrixXd& sol, int i)
{
	VectorXd Fi, solt;
	Fi.resize(5);

	if ((i!=0)&&(i!=_N-1))
	{
		solt = 0.5*( 3*sol.col(i)-sol.col(i+1) );

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( sol.col(i-1)+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
	}
	else if (i==0)
	{
		solt = 0.5*( 3*sol.col(i)-sol.col(i+1) );

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( _Ul+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
	}
	else
	{
		solt = 0.5*( 3*sol.col(i)-_Ur );

		Fi(0) = solt(1);
		Fi(1) = ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) = ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) = solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) = ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));

		solt = 0.5*( sol.col(i-1)+sol.col(i) );

		Fi(0) += solt(1);
		Fi(1) += ( solt(1)*solt(1) - solt(3)*solt(3) )/solt(0) + 0.5*g*solt(0)*solt(0);
		Fi(2) += ( solt(1)*solt(2) - solt(3)*solt(4) )/solt(0);
		Fi(3) += solt(3)*(sol(1,i)/sol(0,i));
		Fi(4) += ( solt(1)*solt(4) - solt(2)*solt(3) )/solt(0) + solt(3)*(sol(2,i)/sol(0,i));
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
	else if(abs(w-v)<5e-13)
		theta = 0;
	else
		theta = (v-u)/(w-u);

	if((theta<0.5)||(theta >2))
		theta = 0;

	phi = (theta + abs(theta))/(1.+theta);

	return phi;
}


// -----------------------------------------
//   Relaxation method
// -----------------------------------------
WRS::WRS(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{}

double WRS::vp_b(const Eigen::MatrixXd& sol, int i)
{

}

void WRS::BuildF(const double& t, const Eigen::MatrixXd& sol)
{

}

VectorXd WRS::Flux_L(const Eigen::MatrixXd& sol, int i)
{

}

VectorXd WRS::Flux_R(const Eigen::MatrixXd& sol, int i)
{

}

// Sauvegarde la solution
void SpaceScheme::SaveSol(const Eigen::MatrixXd& sol, int n)
{
	string name_file = _results + "/solution_" + std::to_string(n) + ".dat";

	assert((sol.cols() == _N) && "The size of the solution matrix is not the same than the number of meshes !");
	assert((sol.rows() == 5) && "The size of the solution matrix is not the same than the number of variables !");

	double dx = 1./_N;

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(14);

	for (size_t i = 0; i < _N; i++)
	{
		solution << i*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																				<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
		solution << (i+1)*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																						<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
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
	solution.precision(7);

	solution << 0. << " " << _Ul(0) << " " << _Ul(1)/_Ul(0) << " " << _Ul(2)/_Ul(0)
																	<< " " << _Ul(3)/_Ul(0) << " " << _Ul(4)/_Ul(0) << endl;
	for (size_t i = 0; i < _N; i++)
		solution << (i+0.5)*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																							<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;

  solution << 1. << " " << _Ur(0) << " " << _Ur(1)/_Ur(0) << " " << _Ur(2)/_Ur(0)
																	<< " " << _Ur(3)/_Ur(0) << " " << _Ur(4)/_Ur(0);

	solution.close();
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
	bruh.precision(7);

	bruh << 0. << " " << _Ul(0) << " " << _Ul(1)/_Ul(0) << " " << _Ul(2)/_Ul(0)
															<< " " << _Ul(3)/_Ul(0) << " " << _Ul(4)/_Ul(0) << endl;
	for (size_t i = 0; i < N; i++)
		bruh << (i+0.5)*dx << " " << exact_sol(0,i) << " " << exact_sol(1,i) << " " << exact_sol(2,i)
																								<< " " << exact_sol(3,i) << " " << exact_sol(4,i) << endl;
  bruh << 1. << " " << _Ur(0) << " " << _Ur(1)/_Ur(0) << " " << _Ur(2)/_Ur(0)
															<< " " << _Ur(3)/_Ur(0) << " " << _Ur(4)/_Ur(0);
	bruh.close();
	// --------------------

	error.setZero(5);; eL2.setZero(5);
	k=0;
	for (int i = 0 ; i < N ; i++)
	{
		if((i+1)*dx <= (k+1)*_dx)
		{
			error(0) += dx*(exact_sol(0,i) - sol(0,k))*(exact_sol(0,i) - sol(0,k));
			for(int j = 1; j<5; j++)
				error(j) += dx*(exact_sol(j,i) - sol(j,k)/sol(0,k))*(exact_sol(j,i) - sol(j,k)/sol(0,k));
		}
		else
		{
			error(0) += ((k+1)*_dx  - i*dx)*(exact_sol(0,i) - sol(0,k))*(exact_sol(0,i) - sol(0,k));
			error(0) += ((i+1)*dx - (k+1)*_dx)*(exact_sol(0,i) - sol(0,k+1))*(exact_sol(0,i) - sol(0,k+1));
			for(int j = 1; j<5; j++)
			{
				error(j) += ((k+1)*_dx  - i*dx)*(exact_sol(j,i) - sol(j,k)/sol(0,k))*(exact_sol(j,i) - sol(j,k)/sol(0,k));
				error(j) += ((i+1)*dx - (k+1)*_dx)*(exact_sol(j,i) - sol(j,k+1)/sol(0,k+1))*(exact_sol(j,i) - sol(j,k+1)/sol(0,k+1));
			}
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
