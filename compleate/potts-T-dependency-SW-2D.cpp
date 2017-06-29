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

int L = 32;
int N = L*L;                          // number of spins in x and y
int n_q = 2;                          // number of state 
// S_i = 0, 1,.., n_q-1

double J_p = 1.0;                   // spin-spin coupling constant
double T;                           // temperature in units of k_B/J_0
vector < vector<int> > S( L, vector<int>(L,0) );  
vector < vector<int> > h( L, vector<int>(L,0) );  
/*SW*/
vector < vector<int> > link_x( L, vector<int>(L,0) );
vector < vector<int> > link_y( L, vector<int>(L,0) );
vector <int> set(N,-1); // Each element remember their parent's index.
vector<int> leaves;
vector <int> level(N,0);	// it lakes 0 or 1.

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

int get_root(int x){
	int root_x;
	if(set[x]<0){ root_x = x; }
	if(set[x]>=0){ root_x = get_root(set[x]); }
	return root_x;
}

void myunion(int root_x, int root_y){
		/* Recall that set[root] remember their depth.*/
		/* Longer root is more robast.*/
		if( (-set[root_x]) > (-set[root_y]) ){
			set[root_y] = root_x;
			level[root_x] = level[root_y]+1 ;
		}
		if( (-set[root_x]) < (-set[root_y]) ){
			set[root_x] = root_y;
			level[root_y] = level[root_x]+1;	
		}	
		if( (-set[root_x]) == (-set[root_y]) ){
			// root_x is chosen as new root.
			set[root_x] -= 1; // Added 1 to the depth of new root.
			set[root_y] = root_x; 
			level[root_x] = level[root_y]+1 ;
		}
}

bool myfind(int x, int y){
	bool flag=false;
	int root_x, root_y;
	root_x = get_root(x);	
	root_y = get_root(y);
	if(root_x != root_y){
		myunion(root_x, root_y);	
		flag = true;
	}
	//print_set_and_level();
	return flag;
}

// input: parent, output: their leaves.
void get_leaves(int i){
	for(int j=0;j<N;++j){
		if(set[j]==i){
			leaves.push_back(j);
			if(level[j]>0){ 
				get_leaves(j);
			}
	}}
}

int heat_prop(int state){
	return rand()%n_q;
}
// this should be completed:w
int nested_prop(int state){
	return 0;
}
	
void SW(double T){
/**** Link Cutting or Freezing Process. *******/	
	//double p_f = 1.0 - exp(-2*J_p/T); // probability of frozen-bond. 
	double p_f = 1.0 - exp(-J_p/T); // probability of frozen-bond. 
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			set[i*L+j]=-1;// initial stete: all nodes are considerd to be root.
			level[i*L+j]= 0;// initial stete: they do not have any leaves.
			link_x[i][j] = 0; link_y[i][j] = 0;	
	}	} // Initialize. This vector tells index of cluster.
	int x_plus,y_plus;
/*******Finding the bond amoung the lattice.*********/
	for(int x=0;x<L;x++){
		for(int y=0;y<L;y++){
			x_plus=(x+1)%L; y_plus=(y+1)%L;
			// This condition can be located in the outside of this function.
			bool condtion_x = kroneckerDelta(S[x][y],S[x_plus][y] )  > 0; // ==1
			bool condtion_y = kroneckerDelta(S[x][y],S[x][y_plus] )  > 0; // ==1
			
			if( ( std_rand()<p_f ) && condtion_x ){ 
			link_x[x][y]=1; 
			myfind( x*L + y, x_plus*L + y );
			}
			if( ( std_rand()<p_f ) && condtion_y ){ 
			link_y[x][y]=1; 
			myfind( x*L + y, x*L + y_plus );
			}
		}}
/**** Update by using union-find algo. *******/	
	for(int k=0;k<N;++k){
		if(set[k]<0){ 
			int state_of_cluster;
			double r=std_rand();
			state_of_cluster = rand()%n_q; // Decited the color of the cluster.
			leaves.clear();	// initialize.	
		get_leaves(k);
		int i1,i2;
		i1 = int( double(k) /double(L) );
		i2 = k % L;
		S[i1][i2] =state_of_cluster;
		
		int size_cluster = leaves.size();
		if( size_cluster > 0 ){
		/* Update-step; Spins that belong to the same cluster take same value.*/	
			for(int l=0;l<size_cluster;++l){
				i1 = int( double(leaves[l]) /double(L) );
				i2 = leaves[l]%L;
				S[i1][i2] =state_of_cluster;
			}
		}}}
}

int main()
{
    cout << " Monte Carlo Simulation of the 2-dimensional Potts Model, using SW-algo.\n"
         << " ----------------------------------------------------\n"
         << " " << L*L << " spins on " << L << "x" << L << " lattice,"<<endl;
    //init_h(); // introduce external field.

    int n_T = 200;
    double T_min=0.1, T_max=3.1;
    double dt = (T_max-T_min)/double(n_T);
    int equilibration_steps = 2000;
    int production_steps = 120000;
    int np = 60;
    cout << "# T " << " <e> " << " <e^2> " << " <e^4> " << " <m> " 
	 << " <m^2> "  << " <m^4> " << " c " << " chi "<< " %accepts "<< endl;
    
    for(int k =0;k<n_T;++k){
	    T = T_min + k * dt; 
    
	    init_S();
	    for (int step = 0; step < equilibration_steps; ++step)
		SW(T);		

	    double e_av = 0, e_sqd_av = 0, e_qua_av = 0;      // average E per spin, squared average
	    double m_av = 0, m_sqd_av = 0, m_qua_av = 0;      // average E per spin, squared average
	    for (int step = 0; step < production_steps; ++step) {
		double num_spins = N;
		if(step%np == 0){
		SW(T);		
		double e_step = E() / num_spins;
		double m_step = M() / num_spins;
		e_av += e_step;
		e_sqd_av += e_step * e_step;
		e_qua_av += e_step * e_step * e_step * e_step;
		m_av += abs(m_step);
		m_sqd_av += m_step * m_step;
		m_qua_av += m_step * m_step * m_step * m_step;
		}
	    }
	    double per_step= double(np) / double(production_steps);
	    e_av *= per_step;
	    e_sqd_av *= per_step;
	    e_qua_av *= per_step;
	    m_av *= per_step;
	    m_sqd_av *= per_step;
	    m_qua_av *= per_step;
	    double c = (e_sqd_av - e_av * e_av) / (T * T);
	    double chi = (m_sqd_av - m_av * m_av) / T;
	     
	    cout << T <<" " << e_av << " " << e_sqd_av << " " << e_qua_av << " " 
		 << m_av << " " << m_sqd_av << " " << m_qua_av << " " 
		 << c << " " << chi << endl; 
       }
    return 0;
}
