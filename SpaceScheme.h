#ifndef _SPACESCHEME_H
#include <string>
#include "Dense"
#include "Sparse"
#include "DataFile.h"

#define g 9.8

class SpaceScheme {
	protected:
		int _N;
		const std::string _results;

		Eigen::MatrixXd _U, _F;
		Eigen::VectorXd _Ul, _Ur;
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

	  // Condition Initiale au centre des triangles
	  Eigen::MatrixXd InitialCondition();

	  // Solution exacte au centre des triangles
	  Eigen::MatrixXd ExactSolution(double t);

	  // Sauvegarde la solution
	  void SaveSol(const Eigen::MatrixXd& sol, int n);

		// Calcul de l'erreur en norme L2
		void ComputeError(const double t);
};

class Rusanov : public SpaceScheme
{
  public:
		// Constructeur
		Rusanov(DataFile* data_file);

    // Une étape du schéma en temps
    void BuildF(const double& t, const Eigen::MatrixXd& sol);
};

class WRS : public SpaceScheme
{
  public:
		// Constructeur
		WRS(DataFile* data_file);

    // Une étape du schéma en temps
    void BuildF(const double& t, const Eigen::MatrixXd& sol);
};


#define _SCPACESCHEME_H
#endif
