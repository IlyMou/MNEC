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

Rusanov::Rusanov(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{}

WRS::WRS(DataFile* data_file) : SpaceScheme::SpaceScheme(data_file)
{}

// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
void Rusanov::BuildF(const double& t, const Eigen::MatrixXd& sol)
{
	double b,bt;
	for(int i=1; i <_N-2; i++)
	{
		_F(0,i) = 0.5*( sol(1,i+1) - sol(1,i-1) );
		_F(1,i) = 0.5*( (sol(1,i+1)*sol(1,i+1)/sol(0,i+1) + 0.5*g*sol(0,i+1)*sol(0,i+1) - sol(3,i+1)*sol(3,i+1)/sol(0,i+1))
									- (sol(1,i-1)*sol(1,i-1)/sol(0,i-1) + 0.5*g*sol(0,i-1)*sol(0,i-1) - sol(3,i-1)*sol(3,i-1)/sol(0,i-1)) ) ;
		_F(2,i) = 0.5*( (sol(1,i+1)*sol(2,i+1)/sol(0,i+1) - sol(3,i+1)*sol(4,i+1)/sol(0,i+1))
		 							-	(sol(1,i-1)*sol(2,i-1)/sol(0,i-1) - sol(3,i-1)*sol(4,i-1)/sol(0,i-1)) );
		_F(3,i) = 0.5*( sol(3,i+1) - sol(3,i-1) )*(sol(1,i)/sol(0,i)) ;
		_F(4,i) = 0.5*( ( sol(1,i+1)*sol(4,i+1) - sol(2,i+1)*sol(3,i+1) )/sol(0,i+1)
		 							- ( sol(1,i-1)*sol(4,i-1) - sol(2,i-1)*sol(3,i-1) )/sol(0,i-1) ) + 0.5*( sol(3,i+1) - sol(3,i-1) )*(sol(2,i)/sol(0,i));

		bt = max(sol(1,i)+sqrt(sol(3,i)*sol(3,i)+g*sol(0,i)),sol(1,i+1)+sqrt(sol(3,i+1)*sol(3,i+1)+g*sol(0,i+1)));
		b  = max(abs(sol(1,i)-sqrt(sol(3,i)*sol(3,i)+g*sol(0,i))),abs(sol(1,i+1)-sqrt(sol(3,i+1)*sol(3,i+1)+g*sol(0,i+1))));
		b  = max(b,bt);
		_F(0,i) -= 0.5*b*(sol(0,i+1)-2*sol(0,i)+sol(0,i-1));
		_F(1,i) -= 0.5*b*(sol(1,i+1)-2*sol(1,i)+sol(1,i-1));
		_F(2,i) -= 0.5*b*(sol(2,i+1)-2*sol(2,i)+sol(2,i-1));
		_F(3,i) -= 0.5*b*(sol(3,i+1)-2*sol(3,i)+sol(3,i-1));
		_F(4,i) -= 0.5*b*(sol(4,i+1)-2*sol(4,i)+sol(4,i-1));
	}
	// Bord gauche
	_F(0,0) = 0.5*( sol(1,1) - _Ul(1) );
	_F(1,0) = 0.5*( (sol(1,1)*sol(1,1)/sol(0,1) + 0.5*g*sol(0,1)*sol(0,1) - sol(3,1)*sol(3,1)/sol(0,1))
								- (_Ul(1)*_Ul(1)/_Ul(0) + 0.5*g*_Ul(0)*_Ul(0) - _Ul(3)*_Ul(3)/_Ul(0)) ) ;
	_F(2,0) = 0.5*( (sol(1,1)*sol(2,1)/sol(0,1) - sol(3,1)*sol(4,1)/sol(0,1))
								-	(_Ul(1)*_Ul(2)/_Ul(0) - _Ul(3)*_Ul(4)/_Ul(0)) );
	_F(3,0) = 0.5*( sol(3,1) - _Ul(3) )*(sol(1,0)/sol(0,0)) ;
	_F(4,0) = 0.5*( ( sol(1,1)*sol(4,1) - sol(2,1)*sol(3,1) )/sol(0,1)
								- ( _Ul(1)*_Ul(4) - _Ul(2)*_Ul(3) )/_Ul(0) ) + 0.5*( sol(3,1) - _Ul(3) )*(sol(2,0)/sol(0,0));

	bt = max(sol(1,0)+sqrt(sol(3,0)*sol(3,0)+g*sol(0,0)),sol(1,1)+sqrt(sol(3,1)*sol(3,1)+g*sol(0,1)));
	b  = max(abs(sol(1,0)-sqrt(sol(3,0)*sol(3,0)+g*sol(0,0))),abs(sol(1,1)-sqrt(sol(3,1)*sol(3,1)+g*sol(0,1))));
	b  = max(b,bt);
	_F(0,0) -= 0.5*b*(sol(0,1)-2*sol(0,0)+_Ul(0));
	_F(1,0) -= 0.5*b*(sol(1,1)-2*sol(1,0)+_Ul(1));
	_F(2,0) -= 0.5*b*(sol(2,1)-2*sol(2,0)+_Ul(2));
	_F(3,0) -= 0.5*b*(sol(3,1)-2*sol(3,0)+_Ul(3));
	_F(4,0) -= 0.5*b*(sol(4,1)-2*sol(4,0)+_Ul(4));

	// Bord droit
	_F(0,_N-1) = 0.5*( _Ur(1) - sol(1,_N-2) );
	_F(1,_N-1) = 0.5*( (_Ur(1)*_Ur(1)/_Ur(0) + 0.5*g*_Ur(0)*_Ur(0) - _Ur(3)*_Ur(3)/_Ur(0))
								- (sol(1,_N-2)*sol(1,_N-2)/sol(0,_N-2) + 0.5*g*sol(0,_N-2)*sol(0,_N-2) - sol(3,_N-2)*sol(3,_N-2)/sol(0,_N-2)) ) ;
	_F(2,_N-1) = 0.5*( (_Ur(1)*_Ur(2)/_Ur(0) - _Ur(3)*_Ur(4)/_Ur(0))
								-	(sol(1,_N-2)*sol(2,_N-2)/sol(0,_N-2) - sol(3,_N-2)*sol(4,_N-2)/sol(0,_N-2)) );
	_F(3,_N-1) = 0.5*( _Ur(3) - sol(3,_N-2) )*(sol(1,_N-1)/sol(0,_N-1)) ;
	_F(4,_N-1) = 0.5*( ( _Ur(1)*_Ur(4) - _Ur(2)*_Ur(3) )/_Ur(0)
								- ( sol(1,_N-2)*sol(4,_N-2) - sol(2,_N-2)*sol(3,_N-2) )/sol(0,_N-2) ) + 0.5*( _Ur(3) - sol(3,_N-2) )*(sol(2,_N-1)/sol(0,_N-1));

	bt = max(sol(1,_N-1)+sqrt(sol(3,_N-1)*sol(3,_N-1)+g*sol(0,_N-1)),_Ur(1)+sqrt(_Ur(3)*_Ur(3)+g*_Ur(0)));
	b  = max(abs(sol(1,_N-1)-sqrt(sol(3,_N-1)*sol(3,_N-1)+g*sol(0,_N-1))),abs(_Ur(1)-sqrt(_Ur(3)*_Ur(3)+g*_Ur(0))));
	b  = max(b,bt);
	_F(0,_N-1) -= 0.5*b*(_Ur(0)-2*sol(0,_N-1)+sol(0,_N-2));
	_F(1,_N-1) -= 0.5*b*(_Ur(1)-2*sol(1,_N-1)+sol(1,_N-2));
	_F(2,_N-1) -= 0.5*b*(_Ur(2)-2*sol(2,_N-1)+sol(2,_N-2));
	_F(3,_N-1) -= 0.5*b*(_Ur(3)-2*sol(3,_N-1)+sol(3,_N-2));
	_F(4,_N-1) -= 0.5*b*(_Ur(4)-2*sol(4,_N-1)+sol(4,_N-2));
}

void WRS::BuildF(const double& t, const Eigen::MatrixXd& sol)
{

}

// Solution exacte au centre des triangles
MatrixXd SpaceScheme::ExactSolution(const double t)
{
	MatrixXd exact_sol(5,_N);
	for (int i = 0 ; i < _N ; i++)
	{
		exact_sol(0,i) = 0;
		exact_sol(1,i) = 0;
		exact_sol(2,i) = 0;
		exact_sol(3,i) = 0;
		exact_sol(4,i) = 0;
	}
	return exact_sol;
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
	solution.precision(7);

	for (size_t i = 0; i < _N; i++)
	{
		solution << i*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																				<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
		solution << (i+1)*dx << " " << sol(0,i) << " " << sol(1,i)/sol(0,i) << " " << sol(2,i)/sol(0,i)
																						<< " " << sol(3,i)/sol(0,i) << " " << sol(4,i)/sol(0,i) << endl;
	}
  solution << endl;
	solution.close();

	if (_norm_l2 ==  "yes")
	{
		if (n%10 == 0)
		cout << "Error L2" << endl;
	}
}

void SpaceScheme::ComputeError(const double t)
{
	double error = 0;
	cout << t*error << endl;
}

#define _SCPACESCHEME_CPP
#endif
