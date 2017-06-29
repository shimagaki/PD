//  Two-dimensional Potts Model
// iAskin , Juius; Teller , Edward(1943). "Statistics of Two-Dimensional Lattices With Four Componets"
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

int L = 16;
int N = L*L;                          // number of spins in x and y
int n_q = 2;                          // number of state 
// S_i = 0, 1,.., n_q-1

double J_p = 1.0;                   // spin-spin coupling constant
double T;                           // temperature in units of k_B/J_0
vector < vector<int> > S( L, vector<int>(L,0) );  
vector < vector<int> > h( L, vector<int>(L,0) );  

void init_S(){
	for(int i=0;i<L;++i){
		for(int j=0;j<L;++j){
		S[i][j] = rand()%n_q;
	}}
}

void init_h(){
	for(int i=0;i<L;++i){
		for(int j=0;j<L;++j){
		h[i][j] =  1;
	}}
}

int kroneckerDelta(int u, int v){
	int w;
	if( u==v )
		w = 1;
	if( u!=v )
		w = 0;
	return w;
}

double E(int x, int y)
{
    int x_plus = x + 1, x_minus = x - 1, y_plus = y + 1, y_minus = y - 1;
    if (x_plus  ==  L)  x_plus  = 0;
    if (x_minus == -1)  x_minus = L - 1;
    if (y_plus  ==  L)  y_plus  = 0;
    if (y_minus == -1)  y_minus = L - 1;

    // contribution to energy from 4 bonds at site (x,y)
    double E = -J_p * (
        kroneckerDelta(S[x][y] , S[x_plus][y]) +
        kroneckerDelta(S[x][y] , S[x_minus][y]) +
	kroneckerDelta(S[x][y] , S[x][y_plus]) +
	kroneckerDelta(S[x][y] , S[x][y_minus]) );
    return E;
}

double E()
{
    // total energy
    double E_sum = 0;
    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y)
            E_sum += E(x, y);
    return E_sum / 2;
}
//  m := q/(q-1) <(S_i -1/q)>, this formulation is correspond to the case of Ising(q=2).
double M(){
	double m=0;
	for(int x=0;x<L;++x){
		for(int y=0;y<L;++y){
		m += (S[x][y] - 1.0/n_q);	
	}}
	return m * n_q/(n_q-1);
}

int heat_prop(int state){
	return rand()%n_q;
}
// this should be completed:w
int nested_prop(int state){
	return 0;
}

bool metropolis_step(           // true/false if step accepted/rejected
    int x, int y,               // change spin at size x,y
    double delta = 1.0 )	// maximum step size in radians 
{
    // find the local energy at site (x,y)
    double E_xy = E(x, y);

    // save the value of spin at (x,y)
    double S_save = S[x][y];

    // trial Metropolis step for spin at site (x,y)
    S[x][y] = heat_prop(S[x][y]);

    double E_xy_trial = E(x, y);

    // Metropolis test
    bool accepted = false;
    double dE = E_xy_trial - E_xy;
    double w = exp(- dE / T);
    if (w > std_rand())
        accepted = true;
    else
        S[x][y] = S_save;
    return accepted;
}

int monte_carlo_sweep()
{
    int acceptances = 0;
    // randomized sweep
    for (int step = 0; step < N; ++step) {
        // choose a site at random
        int x = int(L * std_rand());
        int y = int(L * std_rand());
        if (metropolis_step(x, y))
            ++acceptances;
    }
    return acceptances;
}

int main()
{
    cout << "# Monte Carlo Simulation of the 2-dimensional Potts Model\n"
         << "# ----------------------------------------------------\n"
         << "# " << L*L << " spins on " << L << "x" << L << " lattice,"<<endl;
    //init_h(); // introduce external field.

    int n_T = 200;
    double T_min=0.1, T_max=3.1;
    double dt = (T_max-T_min)/double(n_T);
    int equilibration_steps = 2000;
    int production_steps = 120000;
    int np = 60;
    cout << "# T " << " <e> " << " <e^2> " << " <m> " 
	 << " <m^2> " << " c " << " chi "<< " %accepts "<< endl;
    for(int k =0;k<n_T;++k){
	    T = T_min + k * dt; 
     
	    init_S();
	    for (int step = 0; step < equilibration_steps; ++step)
		monte_carlo_sweep();

	    double e_av = 0, e_sqd_av = 0;      // average E per spin, squared average
	    double m_av = 0, m_sqd_av = 0;      // average E per spin, squared average
	    double accept_percent = 0;
	    for (int step = 0; step < production_steps; ++step) {
		double num_spins = N;
		if(step%np == 0){
		accept_percent += monte_carlo_sweep() / num_spins;
		double e_step = E() / num_spins;
		double m_step = M() / num_spins;
		e_av += e_step;
		e_sqd_av += e_step * e_step;
		m_av += m_step;
		m_sqd_av += m_step * m_step;
		}
	    }
	    double per_step= double(np) / double(production_steps);
	    accept_percent *= per_step * 100.0;
	    e_av *= per_step;
	    e_sqd_av *= per_step;
	    m_av *= per_step;
	    m_sqd_av *= per_step;
	    double c = (e_sqd_av - e_av * e_av) / (T * T);
	    double chi = (m_sqd_av - m_av * m_av) / T;
    
    cout << T <<" " << e_av << " " << e_sqd_av <<" " 
    	 << m_av << " " << m_sqd_av <<" " 
    	 << c << " " << chi <<" " << accept_percent << endl; 
    }  
     return 0;
}
