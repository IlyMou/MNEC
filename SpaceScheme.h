#ifndef _SPACESCHEME_H
#include <string>
#include "Dense"
#include "Sparse"
#include "DataFile.h"

#define g 9.8

class SpaceScheme {
	protected:
		int _N;
		double _bmax, _x0;
		const std::string _results;

		Eigen::MatrixXd _U, _F;
		Eigen::VectorXd _Ul, _Ur, Ul, Ur;
 		std::string _norm_l2;

		DataFile* _data_file;

	public:
		// Constructeur
		SpaceScheme(DataFile* data_file);

		// calcul de V au centre des triangles
	  void BuildVelocity(const double t);

		// Construit le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
	  virtual void BuildF(const double& t, const Eigen::MatrixXd& sol) = 0;

		// Renvoie le vecteur f = F(u,t) (EDO : du/dt = F(u,t))
	  Eigen::MatrixXd& GetF() {return _F;};

		virtual Eigen::VectorXd Flux_R(const Eigen::MatrixXd& sol, int i) = 0;
		virtual Eigen::VectorXd Flux_L(const Eigen::MatrixXd& sol, int i) = 0;

		virtual double vp_b(const Eigen::MatrixXd& sol, int i);

	  // Condition Initiale au centre des triangles
	  Eigen::MatrixXd InitialCondition();

	  // Solution exacte au centre des triangles
	  Eigen::MatrixXd ExactSolution();

	  // Sauvegarde la solution
	  void SaveSol(const Eigen::MatrixXd& sol, int n);

		// Sauvegarde la solution
		void SaveSol(const Eigen::MatrixXd& sol);

		// Calcul de l'erreur en norme L2
		void ComputeError(const Eigen::MatrixXd& sol);

		// Renvoie la valeur maximal de vitesse
		double& Getbmax() {return _bmax;};
};

class Rusanov1 : public SpaceScheme
{
  public:
		// Constructeur
		Rusanov1(DataFile* data_file);

    // Une étape du schéma en temps
    void BuildF(const double& t, const Eigen::MatrixXd& sol);

		Eigen::VectorXd Flux_R(const Eigen::MatrixXd& sol, int i);

		Eigen::VectorXd Flux_L(const Eigen::MatrixXd& sol, int i);
};

class Rusanov2 : public SpaceScheme
{
	private:
		bool _stab;
		SpaceScheme* _rO1;

  public:
		// Constructeur
		Rusanov2(DataFile* data_file);

    // Une étape du schéma en temps
    void BuildF(const double& t, const Eigen::MatrixXd& sol);

		Eigen::VectorXd Flux_R(const Eigen::MatrixXd& sol, int i);

		Eigen::VectorXd Flux_L(const Eigen::MatrixXd& sol, int i);

		double limPente(const Eigen::MatrixXd& sol, int var, int i);
};

class WRS : public SpaceScheme
{
  public:
		// Constructeur
		WRS(DataFile* data_file);

    // Une étape du schéma en temps
    void BuildF(const double& t, const Eigen::MatrixXd& sol);

		Eigen::VectorXd Flux_R(const Eigen::MatrixXd& sol, int i);

		Eigen::VectorXd Flux_L(const Eigen::MatrixXd& sol, int i);

		double vp_b(const Eigen::MatrixXd& sol, int i);
};


#define _SCPACESCHEME_H
#endif
