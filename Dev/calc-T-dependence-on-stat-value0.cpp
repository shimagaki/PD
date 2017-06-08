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

int N = 32;                          // number of spins in x and y
double J_0 = 1.0;                   // spin-spin coupling constant
double f = 0.0;                     // uniform frustration constant
double T;                           // temperature in units of k_B/J_0
int n_vortex_plus;
int n_vortex_minus;

vector< vector<double> > theta;     // spin variable at each lattice site
vector< vector<double> > psi;       //
vector <int> vortex_x; // find vortex location;
vector <int> vortex_y; // find vortex location;
vector <double> T_list(100,0); // find vortex location;
vector <double> correlation(N/2.0,0); // find vortex location;
vector <double> stiffnes_of_vortex(N*N,0); // average of stiffnes under fixed vortex number, which is smaller than total spin..;
vector <double> stiffnes_of_vortex_add_indx(N*N,0); // average of stiffnes under fixed vortex number, which is smaller than total spin..;

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


double M()
{	double mx=0.0, my=0.0,m;
    // total energy
    for (int x = 0; x < N; ++x){
        for (int y = 0; y < N; ++y){
            mx += cos(theta[x][y]);
            my += sin(theta[x][y]);
	}}
   m = sqrt(mx*mx+my*my);
	   return m;
}

//Helicity_Modulas
double HM()
{	double hm=0,hm_c=0,hm_s=0;
    vector<double> hm_s_x(N-1,0);
    vector<double> hm_s_y(N-1,0);
	// does not take into account the variables of boundary. 
    for (int x = 1; x < N-1; ++x){
        for (int y = 1; y < N-1; ++y){
    		int x_plus = x+1, y_plus = y+1;
   		hm_c += cos(theta[x][y]-theta[x_plus][y])+cos(theta[x][y]-theta[x][y_plus]);
		hm_s_x[y-1] += sin(theta[x][y]-theta[x_plus][y]);
		hm_s_y[x-1] += sin(theta[x][y]-theta[x][y_plus]); 
	}
    }
    	for(int l=1;l<N-1;++l){
		hm_s += hm_s_x[l-1]*hm_s_x[l-1];
		hm_s += hm_s_y[l-1]*hm_s_y[l-1];
	}
	hm = hm_c - hm_s/T;	
	hm /=((N-1.0)*(N-1.0));	
	hm /=2.0;	
	return hm;
}

void  calc_correlation_length(){
    // total energy
    for(int x = 0; x < N; ++x){
        for(int y = 0; y < N; ++y){
		for(int l=1;l<correlation.size();l++){ // l is smaller than half of the lengath
			correlation[l]+= 0.25*(cos(theta[x][y]-theta[(x+l)%N][y])
				+cos(theta[x][y]-theta[x][(y+l)%N])	
				+cos(theta[x][y]-theta[(x-l+N)%N][y])	
				+cos(theta[x][y]-theta[x][(y-l+N)%N]) );	
		}
	}}
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
	/*
	double vor,dth12, dth23, dth34, dth45, dth56,dth67,dth78,dth81;
	dth12=abs(acos(cos(theta[i_plus][j_minus]-theta[i_plus][j])));	
	dth23=abs(acos(cos(theta[i_plus][j]-theta[i_plus][j_plus])));
	dth34=abs(acos(cos(theta[i_plus][j_plus]-theta[i][j_plus])));
	dth45=abs(acos(cos(theta[i][j_plus]-theta[i_minus][j_plus])));
	dth56=abs(acos(cos(theta[i_minus][j_plus]-theta[i_minus][j])));
	dth67=abs(acos(cos(theta[i_minus][j]-theta[i_minus][j_minus])));
	dth78=abs(acos(cos(theta[i_minus][j_minus]-theta[i][j_minus])));
	dth81=abs(acos(cos(theta[i][j_minus]-theta[i_plus][j_minus])));
	if(sin(theta[i_plus][j_minus]-theta[i_plus][j])<0){dth12*=-1;}	
	if(sin(theta[i_plus][j]-theta[i_plus][j_plus])<0){dth23*=-1;}	
	if(sin(theta[i_plus][j_plus]-theta[i][j_plus])<0){dth34*=-1;}	
	if(sin(theta[i][j_plus]-theta[i_minus][j_plus])<0){dth45*=-1;}	
	if(sin(theta[i_minus][j_plus]-theta[i_minus][j])<0){dth56*=-1;}	
	if(sin(theta[i_minus][j]-theta[i_minus][j_minus])<0){dth67*=-1;}	
	if(sin(theta[i_minus][j_minus]-theta[i][j_minus])<0){dth78*=-1;}	
	if(sin(theta[i][j_minus]-theta[i_plus][j_minus])<0){dth81*=-1;}	

	vor=dth12+dth23+dth34+dth45+dth56+dth67+dth78+dth81;
	vor = vor/(2.0*M_PI);
	*/	
	
	double vor,dth12, dth23, dth34, dth41;
	dth12 = acos(cos(theta[i][j]-theta[i][j_plus]));
	dth23 = acos(cos(theta[i][j_plus]-theta[i_minus][j_plus]));
	dth34 = acos(cos(theta[i_minus][j_plus]-theta[i_minus][j]));	
	dth41 = acos(cos(theta[i_minus][j]-theta[i][j])); 	
	if(sin(theta[i][j]-theta[i][j_plus])<0){dth12*=-1;}	
	if(sin(theta[i][j_plus]-theta[i_minus][j_plus])<0){dth23*=-1;}	
	if(sin(theta[i_minus][j_plus]-theta[i_minus][j])<0){dth34*=-1;}	
	if(sin(theta[i_minus][j]-theta[i][j])<0){dth41*=-1;}	
	vor = (dth12+dth23+dth34+dth41)/(2.0*M_PI);
	return vor;	
}

// Find locations of vortexes.
void find_vortex(double th){
	int count = 0;
	double temp_th;
	for(int i =0;i<N;++i){
		for(int j=0;j<N;++j){
			//temp_th is -1 to 1
			temp_th = calc_sum_of_neighbor(i,j);
			if(temp_th>th){
			n_vortex_plus +=1;	
			vortex_x.push_back(i);vortex_y.push_back(j);		
			}	
			if(temp_th<-th){
			n_vortex_minus +=1;	
			vortex_x.push_back(i);vortex_y.push_back(j);		
			}}
		}	
	}	

void output_state_to_file(double T, int step, double n_vortex){
	string fname="config-T"+to_string(T)+"-tMC"+to_string(step)+".dat";
	string fname2="vortex-T"+to_string(T)+"-tMC"+to_string(step)+".dat";
	ofstream file(fname);
	ofstream file2(fname2);
	int i,j;
	for(i=0;i<theta.size();i++){
		for(j=0;j<theta[0].size();j++){
			file<<theta[i][j]<<" ";
		}
	file<<endl;
	}
	
	for(int l=0;l<int(n_vortex);l++){
	file2<<vortex_x[l]<<" "<<vortex_y[l]<<endl;
	}
	file.close();
	file2.close();
}

int main()
{
   cout << " Monte Carlo Simulation of the 2-dimensional XY Model\n"
        << " ----------------------------------------------------\n"
        << " " << N*N << " spins on " << N << "x" << N << " lattice,"
        << " uniform frustration f = " << f << endl;

   initialize();
   int equilibration_steps = 4000; // [1].PRB.27.598
   int production_steps = 40002;
   int np = 10; 
   double num_spins = double(N) * double(N);
   double T_min=0.001;
   double dt = 2.0/T_list.size();
   string fname_cor="correlation.dat";
   ofstream file_cor(fname_cor);
   cout<<"T,<e>,<e^2>,c,<m>,<m>^2,chi,<#vortex_plus>,<#vortex_minus>,Accept"<<endl;
   for(int q=0;q<T_list.size();++q){
    T = T_list[q]+T_min;
    //T = 0.4;
    double e_step,s_step;
    for(int k=0;k<T_list.size();++k){T_list[k] = k*dt;} 
    //cout << " done" << endl;
    //cout << " performing " << equilibration_steps
    //     << " equilibration steps ..." << flush;
    for (int step = 0; step < equilibration_steps; ++step){
        monte_carlo_sweep();
    }
    //cout << " done" << endl;
    
    double e_av = 0, e_sqd_av = 0;      // average E per spin, squared average
    double s_av = 0, s_sqd_av = 0;      // average of spin, squared average
    double vortex_step = 0;      // average vortex. 
    double hm_av=0,hm_step=0;
    double accept_percent = 0;
    double vortex_av_plus=0,vortex_av_minus=0;
    for(int l=0;l<correlation.size();++l){correlation[l]=0.0;}
    //cout << " performing " << production_steps
    //     << " production steps ..." << flush;
    for (int step = 0; step < production_steps; ++step) {
    	monte_carlo_sweep();
	if(step%np == 0){ 
        accept_percent += monte_carlo_sweep() / num_spins;
        e_step = E() / num_spins;
        e_av += e_step;
        e_sqd_av += e_step * e_step;
	s_step=M()  / num_spins;
	s_av += s_step; 
        s_sqd_av += s_step * s_step;
    	hm_step= HM();
    	hm_av += hm_step;
	
	vortex_x.clear(); vortex_y.clear(); n_vortex_plus=0; n_vortex_minus=0;	
	find_vortex(0.9);
	vortex_av_plus+=n_vortex_plus;
	vortex_av_minus+=n_vortex_minus;
	stiffnes_of_vortex[int(n_vortex_plus)]+=hm_step;	
	stiffnes_of_vortex_add_indx[int(n_vortex_plus)]+=1;	
	
	calc_correlation_length();
	if(step<10 || step%10000==0){
		output_state_to_file(T,step,vortex_av_plus+vortex_av_minus);
	}
	}
    }
    //cout << " done" << endl;
    double per_step= np / double(production_steps);
    accept_percent *= per_step * 100.0;
    e_av *= per_step;
    e_sqd_av *= per_step;
    s_av *= per_step;
    s_sqd_av *= per_step;
    vortex_av_plus *= per_step;
    vortex_av_minus *= per_step;
    hm_av *= per_step;
    double c = (e_sqd_av - e_av * e_av) / T;
    double chi = (s_sqd_av - s_av * s_av) / T;
    cout<<T<<" "
	<<e_av<<" "
	<<e_sqd_av<<" "
	<<c<<" "
	<<s_av<<" "
	<<s_sqd_av<<" "
	<<chi<<" "
	<<vortex_av_plus<<" "
	<<vortex_av_minus<<" "
	<<hm_av<<" "
	<<accept_percent<<endl;
    
	file_cor<<T<<" ";
	for(int l=0;l<correlation.size();++l){
	file_cor<<correlation[l]*per_step/(double(N*N))<<" ";
	}file_cor<<endl;
    /*
    cout << " T = " << T <<endl;
    cout << " <e> = " << e_av <<endl;
    cout << " <e^2> = " << e_sqd_av << endl;
    cout << " Heat capacity = " << c << endl;
    cout << " <m> = " << s_av <<endl;
    cout << " <m^2> = " << s_sqd_av << endl;
    cout << " Magnetic susceptibility = " << chi << endl;
    cout << " #vortex = " << vortex_av << endl;
    cout << " accepts = " << accept_percent << endl;
    */
   }//end of T loop.
  file_cor.close();
 
    string fname_HM="Helicity_Modulas-of-vortex.dat";
    ofstream file_HM(fname_HM);
    for(int l=0;l<stiffnes_of_vortex.size();++l){
	 if(stiffnes_of_vortex_add_indx[l]>0){
		file_HM<<l<<" "<<l<<" "<<(stiffnes_of_vortex[l]/double(stiffnes_of_vortex_add_indx[l]))<<" "<<endl; 
	}
    }
   	file_HM.close(); 
       return 0;
}  
