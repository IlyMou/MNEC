#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // ----------------------- Fichier de données --------------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();
  // ---------------------------------------------------------------------------

  // ------------------Définition du nombre d'itérations------------------------
  int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())/data_file->Get_dt()));
  data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
  // ---------------------------------------------------------------------------

  // ---------------------------- Résolution  ----------------------------------
  SpaceScheme* pb = NULL;
  if (data_file->Get_numerical_flux_choice() == "wrs")
    pb = new WRS(data_file);
  else if (data_file->Get_order() == 1)
    pb = new Rusanov1(data_file);
  else
    pb = new Rusanov2(data_file);

  TimeScheme* time_scheme = NULL;
  if (data_file->Get_scheme() == "SSPRK2")
    time_scheme = new SSPRK2();
  else if (data_file->Get_scheme() == "RK2")
    time_scheme = new RK2();
  else if (data_file->Get_scheme() == "RK4")
    time_scheme = new RK4();
  else
    time_scheme = new EulerScheme();

  cout << "-------------------------------------------------" << endl;
  cout << "Search h, u and b such that : " << endl;
  cout << "dt h + div hu = 0" << endl;
  cout << "dt hu + div (hu x u - hb x b) + grad (gh^2)/2= 0" << endl;
  cout << "dt hb + div (hb x u - hu x b) + u div hb= 0" << endl;
  cout << "-------------------------------------------------" << endl;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();
  time_scheme->Initialize(data_file, pb); // Initialisation
  time_scheme->SaveSolution(0); // Sauvegarde condition initiale

  int nsave = nb_iterations/4;
  for (int n = 1; n <= nb_iterations; n++) // Boucle en temps
  {
    cout.flush();
    cout << "Progression : " << (double)n/((double)(nb_iterations))*100 << "% \r";

    time_scheme->Advance();
    if((!(n%nsave))||(n==nb_iterations))
      time_scheme->SaveSolution(n);
  }
  time_scheme->SaveSolution();

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  int timeC = floor(t*0.001); bool hmin = false;
  cout << "Time to compute : ";
  if(timeC/3600>=1){
    cout << floor(timeC/3600) << " h ";
    timeC = timeC - floor(timeC/3600)*3600;
    hmin = true;
  }
  if(timeC/60>=1){
    cout << floor(timeC/60) << " min ";
    timeC = timeC - floor(timeC/60)*60;
    hmin = true;
  }
  if(hmin) cout << "and ";
  cout << floor(timeC)+1 << " seconds" << endl;
  // ---------------------------------------------------------------------------


  //------------------- Si on connait la solution exacte -----------------------
  if ((data_file->Get_initial_condition_choice() == "riemann")&&(data_file->Get_norm_l2()=="yes"))
  {
    time_scheme->ComputeError();
  }
  // ---------------------------------------------------------------------------

  delete time_scheme;
  delete pb;
  delete data_file;

  return 0;
}
