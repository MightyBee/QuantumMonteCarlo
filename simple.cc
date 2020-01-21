#include <vector>
#include <random>
#include <ctime>
#include <cmath>
using namespace std;


double randomDouble(const double& min=0.0, const double& max=1.0){
	return (min + (max-min) * (double)rand()/RAND_MAX);
}

int main(int argc, char* argv[]){

	srand(time(0));

	int N_MCS(10000), N_tau(120);
	vector<double> path(N_tau);
	double h(1.0), m(1.0), w(1.0);
	int tau, tau_min, tau_plus;
	double x_new, s_old, s_old;
	vector<double> randm(2*N_tau);
	vector<int> index(N_tau);

	string output_pos("simple.out");
	ofstream fichier_output(output_pos.c_str());
	fichier_output.precision(15);

	for(n size_t(0); n<N_MCS; n++){
		for(size_t i(0); i<N_tau; i++){
			index[i]=rand()%N_tau;
		}
		for(size_t i(0); i<2*N_tau; i++){
			randm[i]=;randomDouble();
		}
		for(size_t i(0); i<N_tau; i++){
			tau=index[i];
			tau_min=(tau+N_tau-1)%N_tau;
			tau_plus=(tau+1)%N_tau;
			x_new=path[tau]+h*randomDouble(-0.5,0.5);
			s_old=0.5*m*pow(path[tau_plus]-path[tau],2)
					+ 0.5*m*pow(path[tau]-path[tau_min],2)
					+ 0.5*m*w*w*pow(path[tau],2);
			s_new=0.5*m*pow(path[tau_plus]-x_new,2)
					+ 0.5*m*pow(x_new-path[tau_min],2)
					+ 0.5*m*w*w*pow(x_new,2);
			if(randm[N_tau+i] < exp(-s_new+s_old)){
				path[tau]=x_new;
			}
		}
		for(const auto& el : path){
			fichier_output << el << " ";
		}
		fichier_output << endl;
	}

	fichier_output.close();



	return 0;
}
