#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name),  _if_N_mesh(false),
_if_t0(false), _if_tfinal(false), _if_dt(false),
_if_scheme(false),_if_initial_condition_choice(false),
 _if_results(false), _if_numerical_flux_choice(false)
{}

void DataFile::ReadDataFile()
{
  double bm,cfl;
  ifstream data_file(_file_name.data());
  if (!data_file.is_open())
  {
    cout << "Unable to open file " << _file_name << endl;
    abort();
  }
  else
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Reading data file " << _file_name << endl;
  }

  string file_line;

  while (!data_file.eof())
  {
    getline(data_file, file_line);
    if (file_line.find("N_mesh") != std::string::npos)
    {
      data_file >> _N_mesh; _if_N_mesh = true;
    }

    if (file_line.find("t0") != std::string::npos)
    {
      data_file >> _t0; _if_t0 = true;
    }

    if (file_line.find("tfinal") != std::string::npos)
    {
      data_file >> _tfinal; _if_tfinal = true;
    }

    if (file_line.find("dt") != std::string::npos)
    {
      data_file >> _dt; _if_dt = true;
    }

    if (file_line.find("initial_condition") != std::string::npos)
    {
      data_file >> _initial_condition_choice; _if_initial_condition_choice = true;
      if (_initial_condition_choice == "gaussian")
      {
        data_file >> _x0 >> _y0 >> _a;
      }
      else if (_initial_condition_choice == "riemann")
      {
        data_file >> _x0 >> _a >> _b >> _c >> _d >> _e >> _f >> _g >> _h >> _i >> _j;

        bm = max(_b-sqrt(_d*_d+g*_a),_g+sqrt(_i*_i+g*_f));
        cfl = bm*_dt*_N_mesh;
      }
      else if (_initial_condition_choice == "rectangular")
      {
        data_file >> _x0 >> _y0 >> _b;
      }
      else
      {
        cout << "Only Riemann problem in 1D and gaussian or rectangular initial conditions in 2D are implemented." << endl;
        abort();
      }
    }

    if (file_line.find("numerical_flux") != std::string::npos)
    {
      data_file >> _numerical_flux_choice; _if_numerical_flux_choice = true;
      if ((_numerical_flux_choice != "rusanov") && (_numerical_flux_choice != "wrs"))
      {
        cout << "Only Rusanov and the 5-Wave Relaxation solver schemes are implemented." << endl;
        abort();
      }
    }

    if (file_line.find("scheme") != std::string::npos)
    {
      data_file >> _scheme; _if_scheme = true;
      if ((_scheme != "ExplicitEuler") && (_scheme != "RK2") && (_scheme != "SSPRK2") && (_scheme != "RK4"))
      {
        cout << "Only Explicit Euler, SSP-RK2, RK2 and RK4 are implemented." << endl;
        abort();
      }
    }

    if (file_line.find("ordre") != std::string::npos)
    {
      data_file >> _order; _if_order = true;
    }

    if (file_line.find("results") != std::string::npos)
    {
      data_file >> _results; _if_results = true;
    }

    if (file_line.find("norm_l2") != std::string::npos)
    {
      data_file >> _norm_l2; _if_norm_l2 = true;
    }
  }

  if (!_if_t0)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.) is used for t0." << endl;
    _t0 = 0.;
  }
  if (!_if_tfinal)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (1) is used for tfinal." << endl;
    _tfinal = 1;
  }
  if (!_if_dt)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (0.001) is used for dt." << endl;
    _dt = 0.001;
  }
  if (!_if_scheme)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default scheme (SSP - Runge Kutta 2) is used." << endl;
    _scheme = "SSPRK2";
  }
  if (!_if_results)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default results folder name (results) is used." << endl;
    _results = "results";
  }
  if (!_if_N_mesh)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default discretization N=100 is used." << endl;
    _N_mesh = 100;
  }
  if (!_if_initial_condition_choice)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default Riemann problem is used with x0=0.5 u_l=1 and u_r=0." << endl;
    _initial_condition_choice = "riemann"; _x0 = 0.5; _a=1; _b=0;
  }
  if (!_if_numerical_flux_choice)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default value (rusanov) is used for the numerical flow." << endl;
    _numerical_flux_choice = "rusanov";
  }
  if (!_if_order)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default order 1 is used." << endl;
    _order = 1;
  }
  if (!_if_norm_l2)
  {
    cout << "-------------------------------------------------" << endl;
    cout << "Beware - The default choice (no) for L2 norm is used." << endl;
    _norm_l2 = "no";
  }
  cout << "-------------------------------------------------" << endl;

  cout << " -- Solver : " << _numerical_flux_choice;
  if(_numerical_flux_choice=="rusanov") cout << "  order " << _order;
  cout << endl;
  cout << " -- Time Scheme : " << _scheme << endl;
  cout << " -- Time step : " << _dt << endl;
  cout << " -- Mesh step : " << _N_mesh << endl;
  cout << " -- CFL : " << cfl << endl;

  if ((((_scheme != "RK2")||(_scheme != "SSPRK2"))&&(cfl>0.25))
  ||((_scheme != "RK4")&&(cfl>0.125))||((_scheme != "ExplicitEuler")&&(cfl>0.5)))
  {
    cout << " ********* Time step too large ********* " << endl << endl;
    abort();
  }

  cout << " -- Initial Condition : " << _initial_condition_choice << endl;
  if(_initial_condition_choice == "riemann")
  {
    cout << "    h = " << _a << "      " << _f << endl;
    cout << "    u = " << _b << "      " << _g << endl;
    cout << "    v = " << _c << "      " << _h << endl;
    cout << "    a = " << _d << "      " << _i << endl;
    cout << "    b = " << _e << "      " << _j << endl;
  }

}

#define _DATA_FILE_CPP
#endif
