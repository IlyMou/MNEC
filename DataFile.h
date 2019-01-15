#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>

#define g 9.8

// Définition de la classe
class DataFile {
private:
  std::string _file_name;

  int _N_mesh, _order, _N_solE;
  double _t0, _tfinal, _dt;
  double _x0, _y0, _a, _b, _c, _d, _e, _f, _g, _h, _i, _j;

  std::string _scheme;
  std::string _initial_condition_choice;
  std::string _numerical_flux_choice;
  std::string _results;
  std::string _norm_l2, _file_solEx;   //calcul de la norme L2 de u

  bool _if_N_mesh;
  bool _if_t0;
  bool _if_tfinal;
  bool _if_dt;
  bool _if_order;
  bool _if_scheme;
  bool _if_initial_condition_choice;
  bool _if_numerical_flux_choice;
  bool _if_results;
  bool _if_norm_l2;

public: // Méthodes et opérateurs de la classe
  DataFile(std::string file_name);
  void ReadDataFile();
  void Adapt_dt(double dt){_dt = dt;};
  int Get_N_mesh() const {return _N_mesh;};
  int Get_order() const {return _order;};
  int Get_N_solE() const {return _N_solE;};
  double Get_t0() const {return _t0;};
  double Get_tfinal() const {return _tfinal;};
  double Get_dt() const {return _dt;};
  double Get_param_x0() const { return _x0;};
  double Get_param_y0() const { return _y0;};
  double Get_param_a() const { return _a;};
  double Get_param_b() const { return _b;};
  double Get_param_c() const { return _c;};
  double Get_param_d() const { return _d;};
  double Get_param_e() const { return _e;};
  double Get_param_f() const { return _f;};
  double Get_param_g() const { return _g;};
  double Get_param_h() const { return _h;};
  double Get_param_i() const { return _i;};
  double Get_param_j() const { return _j;};
  std::string Get_scheme() const {return _scheme;};
  std::string Get_initial_condition_choice() const {return _initial_condition_choice;};
  std::string Get_numerical_flux_choice() const {return _numerical_flux_choice;};
  std::string Get_results() const {return _results;};
  std::string Get_norm_l2() const {return _norm_l2;};
  std::string Get_file_solEx() const {return _file_solEx;};
};

#define _DATA_FILE_H
#endif
