#ifndef _TIME_SCHEME_H

#include "SpaceScheme.h"

class TimeScheme
{
  protected:
    // Nombre de maille
    double _N;
    // Pas de temps
    double _dt;
    // Temps en cours
    double _t;
    // Vecteur initial et vecteur solution
    Eigen::MatrixXd _sol0, _sol;
    // Pointeur vers la classe SpaceScheme
    SpaceScheme* _pb;

  public:
    // Constructeur par défaut
    TimeScheme();
    // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
    // alors le compilateur en génère un implicitement.
    virtual ~TimeScheme();
    // Initialisation de vos différentes variables
    void Initialize(DataFile* data_file, SpaceScheme* pb);
    // Enregistre la solution un fichier
    void SaveSolution(int n);
    // Une étape du schéma en temps
    virtual void Advance() = 0;
    // Permet de récupérer _sol
    const Eigen::MatrixXd & GetIterateSolution() const;
};

class EulerScheme : public TimeScheme
{
  public:
    // Une étape du schéma en temps
    void Advance();
};

class RungeKuttaScheme : public TimeScheme
{
  public:
    // Une étape du schéma en temps
    void Advance();
};

#define _TIME_SCHEME_H
#endif
