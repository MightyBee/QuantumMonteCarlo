#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
using namespace std;

//CIAOOO
//LOL

double randomDouble(const double& min, const double& max); // function ton generate a random (uniform) double between 'min' and 'max'
//double QLagrangian(const vector<vector<double>>& pos, const double& d_tau);
//double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, double const& new_pos);

int main(){
  unsigned int N_part(1), N_slices(20); // number of particles in the system and of (imaginary) time slices respectively
  double tau_final(100.0); // tau<->beta=1/(T*k_B) the so called inverse temperature
  double d_tau(tau_final/N_slices); // imaginary time step
  double mass(1.0), frequency(0.1); // effective mass and angular freqeuncy of the harmonic potential respectively
  double m(mass*d_tau), omega(frequency*d_tau); // dimensionless mass and dimensionless angular frequency respectively
  double pos_min(-10.0), pos_max(+10.0); // initial mininum/maximal position
  double h(10.0); // initial maximum displacement of a point in the path
  vector<vector<double>> positions(N_slices,vector<double>(N_part,0)); // table : each row corresponds to a time slice and contains the coordinates (1D for now) of the N_part particles, or each column corresponds to the path of a particle
  // switch columns and rows, TODO or not TODO ?

  for(auto& Im_t : positions){ // initialize random paths for each particles
    for(auto& pos : Im_t){
      pos=randomDouble(pos_min,pos_max);
    }
  }

  unsigned int MCS(2000); // number of Monte Carlo iterations
  size_t n_stride(1); // output is written in the output file every n_stride iterations
  unsigned int mm(0), mm_plu(0), mm_min(0); // mm : time slice randomly selected during each iteration, mm_plu=mm+1, mm_min=mm-1;
  unsigned int nn(0); // particle randomly selected during each iteration
  double new_r(0.0); // new position proposed
  double s_old(0.0), s_new(0.0); // part of the action that is changed with the new position proposed (old and new values respectively)

<<<<<<< HEAD
  string output("PNGU.out"); // output file
  ofstream fichier_potentiel(output.c_str());
  fichier_potentiel.precision(15);
=======
  string output("caca.out"); // output file
  ofstream fichier_output(output.c_str());
  fichier_output.precision(15);
>>>>>>> cd637479e0bf2f4b1b1add93918324aa8998fceb

  for(size_t i(0); i<MCS; i++){
    for(size_t j(0); j<N_slices; j++){
      mm=rand()%N_slices; // random integer between 0 and N_slices-1
      mm_plu=(mm+N_slices-1)%N_slices; // mm-1 with periodic boundary condition
      mm_min=(mm+1)%N_slices; // mm+1 with periodic boundary condition
      nn=rand()%N_part; // random integer between 0 and N_part-1
      new_r=randomDouble(-1,1)*h; // proposed new position for the particle nn at time slice mm*d_tau

      // as we take the difference of new and old action S_new-S_old, we can consider only the part of the action that is affected by the proposed new position
      s_old=0.5*pow(positions[mm_plu][nn]-positions[mm][nn],2)
           +0.5*pow(positions[mm][nn]-positions[mm_min][nn],2)
           +0.5*pow(positions[mm][nn],2); // harmonic potential
      s_new=0.5*pow(positions[mm_plu][nn]-new_r,2)
           +0.5*pow(new_r-positions[mm_min][nn],2)
           +0.5*pow(new_r,2);

      if(randomDouble(0,1)<exp(s_old-s_new)){ // metropolis acceptance
        positions[mm][nn]=new_r;
      }
    }

    // writing in the output file
    if(!(i%n_stride)){
      for(auto& Im_t : positions){
        for(auto& pos : Im_t){
          fichier_potentiel << pos << " ";
        }
      }
      fichier_potentiel << endl;
    }

  }

  fichier_potentiel.close();

	return 0;
}

double randomDouble(const double& min, const double& max){
  return min + (double)rand()/RAND_MAX*(max-min);
}


/*
double Lagrangian(const vector<vector<double>>& pos, const double& d_tau){
  unsigned int i(0);
  return pow((pos[i+1][0]-pos[i][0])/d_tau,2)+0.5*pow(pos[i][0],2);
}
*/
/*
double diff_QLagrangian(const vector<vector<double>>& pos, const double& d_tau, const unsigned int& m, const unsigned int& n, const double& new_r){
    double Lagr_old(0.0);
    double Lagr_new(0.0);
    vector<vector<double>> pos_new(pos);
    pos_new[m][n]=new_r;
    for(size_t i(0); i<pos.size()-1; i++){
      Lagr_old+=pow((pos[i+1][0]-pos[i][0])/d_tau,2)+0.5*100*pow(pos[i][0],2);
      Lagr_new+=pow((pos_new[i+1][0]-pos_new[i][0])/d_tau,2)+0.5*100*pow(pos_new[i][0],2);
    }
    return Lagr_old-Lagr_new;
}
*/
