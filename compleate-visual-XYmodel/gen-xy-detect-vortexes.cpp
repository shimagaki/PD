//  Two-dimensional XY (planar spin) Model
//  J. Tobochnik and G.V. Chester, Phys. Rev. B 20, 3761 (1979)
//  S. Teitel and C. Jayaprakash, Phys. Rev. B 27, 598-601 (1983)
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

int N = 64;                          // number of spins in x and y
double J_0 = 1.0;                   // spin-spin coupling constant
double f = 0.0;                     // uniform frustration constant
double T;                           // temperature in units of k_B/J_0

vector< vector<double> > theta;     // spin variable at each lattice site
vector< vector<double> > psi;       //
vector <int> vortex_x; // find vortex location;
vector <int> vortex_y; // find vortex location;

void initialize()
{
    theta.resize(N);
    psi.resize(2 * N);
    for (int i = 0; i < N; ++i)
        theta[i].resize(N);
    for (int i = 0; i < 2 * N; ++i)
        psi[i].resize(2 * N);
    if (f != 0.0) {
        // initialize psi
    }
}

double E(int x, int y)
{
    // energy of bonds connected to site (x,y)
    // find the nearest neighbor sites with periodic boundary conditions
    int x_plus = x + 1, x_minus = x - 1, y_plus = y + 1, y_minus = y - 1;
    if (x_plus  ==  N)  x_plus  = 0;
    if (x_minus == -1)  x_minus = N - 1;
    if (y_plus  ==  N)  y_plus  = 0;
    if (y_minus == -1)  y_minus = N - 1;

    // contribution to energy from 4 bonds at site (x,y)
    double E = -J_0 * (
        cos(theta[x][y] - theta[x_plus][y])  +
        cos(theta[x][y] - theta[x_minus][y]) +
        cos(theta[x][y] - theta[x][y_plus])  +
        cos(theta[x][y] - theta[x][y_minus]) );

    return E;
}

double E()
{
    // total energy
    double E_sum = 0;
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
            E_sum += E(x, y);
    return E_sum / 2;
}

bool metropolis_step(               // true/false if step accepted/rejected
    int x, int y,                   // change spin at size x,y
    double delta = 1.0)             // maximum step size in radians
{
    // find the local energy at site (x,y)
    double E_xy = E(x, y);

    // save the value of spin at (x,y)
    double theta_save = theta[x][y];

    // trial Metropolis step for spin at site (x,y)
    theta[x][y] += delta * (2 * std_rand() - 1);

    double E_xy_trial = E(x, y);

    // Metropolis test
    bool accepted = false;
    double dE = E_xy_trial - E_xy;
    double w = exp(- dE / T);
    if (w > std_rand())
        accepted = true;
    else
        theta[x][y] = theta_save;

    return accepted;
}

int monte_carlo_sweep()
{
    int acceptances = 0;
    // randomized sweep
    for (int step = 0; step < N * N; ++step) {
        // choose a site at random
        int x = int(N * std_rand());
        int y = int(N * std_rand());
        if (metropolis_step(x, y))
            ++acceptances;
    }
    return acceptances;
}

double calc_sum_of_neighbor(int i, int j){
	int i_plus = (i+1)%N, j_plus=(j+1)%N
      		,i_minus = (i-1+N)%N, j_minus=(j-1+N)%N;	
	double vor;
	vor =
	/*	
		abs(acos(cos(theta[i_plus][j_minus]-theta[i_plus][j])))	
		+abs(acos(cos(theta[i_plus][j]-theta[i_plus][j_plus])))
		+abs(acos(cos(theta[i_plus][j_plus]-theta[i][j_plus])))
		+abs(acos(cos(theta[i][j_plus]-theta[i_minus][j_plus])))
		+abs(acos(cos(theta[i_minus][j_plus]-theta[i_minus][j])))
		+abs(acos(cos(theta[i_minus][j]-theta[i_minus][j_minus])))
		+abs(acos(cos(theta[i_minus][j_minus]-theta[i][j_minus])))
		+abs(acos(cos(theta[i][j_minus]-theta[i_plus][j_minus])));
	*/
	//	This method is beter than above.	
		+abs(acos(cos(theta[i][j]-theta[i][j_plus])))
		+abs(acos(cos(theta[i][j_plus]-theta[i_minus][j_minus])))
		+abs(acos(cos(theta[i_minus][j_minus]-theta[i_minus][j])))
		+abs(acos(cos(theta[i_minus][j]-theta[i][j])));
	vor = vor/(2.0*M_PI);
	return vor;	
}

// Find locations of vortexes.
int find_vortex(double th){
	int count = 0;
	double temp_th;
	for(int i =0;i<N;++i){
		for(int j=0;j<N;++j){
			temp_th = calc_sum_of_neighbor(i,j);
			if(temp_th>th){
				vortex_x.push_back(i);
				vortex_y.push_back(j);		
				count +=1;
			}		
		}	
	}	
	return count;	
}

void init_vortex(int n){
	double PI2=2.0*M_PI;	
	int M=N-1;	
	int dx=(int)(double(M)/n),dy=(int)(double(M)/n);
	for(int i=1;i<=n;++i){
		int x=dx*i,y=dy*i,sign=pow((-1),n);
		theta[(x+1)%M][(y+1)%M]=PI2*(1.0/8)*sign;	
		theta[(x+1)%M][y]=PI2*(2.0/8)*sign;	
		theta[(x+1)%M][(y+M-1)%M]=PI2*(3.0/8)*sign;	
		theta[x][(y+M-1)%M]=PI2*(4.0/8)*sign;	
		theta[(x+M-1)%M][(y+M-1)%M]=PI2*(5.0/8)*sign;	
		theta[(x+M-1)%M][y]=PI2*(6.0/8)*sign;	
		theta[(x+M-1)%M][(y+1)%M]=PI2*(7.0/8)*sign;	
		theta[x][(y+1)%M]=PI2*(8.0/8)*sign;	
	}
}

int main()
{
    cout << " Monte Carlo Simulation of the 2-dimensional XY Model\n"
         << " ----------------------------------------------------\n"
         << " " << N*N << " spins on " << N << "x" << N << " lattice,"
         << " uniform frustration f = " << f << endl;

    initialize();
    int equilibration_steps = 2000;
    int production_steps = 2002;
    int np = 1000; 
    int step_max=0; 
    T = 10; 
    for (int step = 0; step < 100; ++step){
        monte_carlo_sweep();//set random configuration.
}
    cout << " done" << endl;
    cout << " performing " << equilibration_steps
         << " equilibration steps ..." << flush;
    T = 0.4;
    for (int step = 0; step < equilibration_steps; ++step)
        monte_carlo_sweep();
    cout << " done" << endl;

    double e_av = 0, e_sqd_av = 0;      // average E per spin, squared average
    double accept_percent = 0;
    cout << " performing " << production_steps
         << " production steps ..." << flush;
    for (int step = 0; step < production_steps; ++step) {
	double num_spins = double(N) * double(N);
        accept_percent += monte_carlo_sweep() / num_spins;
        double e_step = E() / num_spins;
        e_av += e_step;
        e_sqd_av += e_step * e_step;
    
    	if(step%np == 0){
		if(step_max<step){step_max=step;}
		string fname="all-spins-tMC"+to_string(step)+".dat";
		string fname2="vortex-tMC"+to_string(step)+".dat";
		ofstream file(fname);
		ofstream file2(fname2);
		int i,j;
		for(i=0;i<theta.size();i++){
			for(j=0;j<theta[0].size();j++){
				file<<theta[i][j]<<" ";
			}
		file<<endl;
		}
		
		vortex_x.clear();vortex_y.clear();	
		int n_vortex=0;
		n_vortex=find_vortex(0.9);
		cout<<"n_vortex="<<n_vortex<<endl;	
		for(int l=0;l<n_vortex;l++){
		file2<<vortex_x[l]<<" "<<vortex_y[l]<<endl;
		//cout<<vortex_x[l]<<","<<vortex_y[l]<<endl;
		}
		file.close();
		file2.close();
    	}
    }
    cout << " done" << endl;

    double per_step= 1.0 / double(production_steps);
    accept_percent *= per_step * 100.0;
    e_av *= per_step;
    e_sqd_av *= per_step;
    cout << " T = " << T << '\t'
         << " <e> = " << e_av << '\t'
         << " <e^2> = " << e_sqd_av << '\t'
         << " %accepts = " << accept_percent << endl;
    double c = (e_sqd_av - e_av * e_av) / (T * T);
    cout << " Heat capacity = " << c << endl;
    
	return 0;
}

