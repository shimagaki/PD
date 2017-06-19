#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

double J_0 = 1.0;	// average strength of the J
double T = 1.0;
vector< vector<double> > J_x;
vector< vector<double> > J_y;

int L = 16;
int N = L*L;
vector <int> set(N,-1); // Each element remember their parent's index.
//	Root rememter their "(-1) x depth", remember every root have negative.
vector < vector<int> > link_x( L, vector<int>(L,0) );
vector < vector<int> > link_y( L, vector<int>(L,0) );
vector<int> leaves;
vector<int> level(N,0);	// it lakes 0 or 1.
vector < vector<int> > S;


inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

void initialize(){
    S.resize(L);
    J_x.resize(L);
    J_y.resize(L);
    for (int i = 0; i < L; ++i){
        S[i].resize(L);
        J_x[i].resize(L);
        J_y[i].resize(L);
	}
}

void init_J_pure_ferro(){
	for(int i=0;i<L;++i){
		for(int j=0;j<L;++j){
		J_x[i][j]= 1.0;
		J_y[i][j]= 1.0;
		
	}}
}

void init_S_neel(){
	for(int i=0;i<L;++i){
		for(int j=0;j<L;++j){
			S[i][j]= pow(-1,i+j);
	}}
}

int get_root(int x){
	cout<<"get_root"<<endl;
	int root_x;
	if(set[x]<0){ root_x = x; }
	if(set[x]>=0){ root_x = get_root(set[x]); }
	return root_x;
}

void myunion(int root_x, int root_y){
		cout<<"myunion"<<endl;
		/* Recall that set[root] remember their depth.*/
		/* Longer root is more robast.*/
		if( (-set[root_x]) > (-set[root_y]) ){
			set[root_y] = root_x;
			level[root_x] = 1 ;
		}
		if( (-set[root_x]) < (-set[root_y]) ){
			set[root_x] = root_y;	
			level[root_y] = 1 ;
		}	
		if( (-set[root_x]) == (-set[root_y]) ){
			set[root_x] -= 1; // Added 1 to the depth of new root.
			set[root_y] = root_x;
			level[root_x] = 1 ;
		}
}

bool myfind(int x, int y){
	cout<<"myfind"<<endl;
	bool flag=false;
	int root_x, root_y;
	root_x = get_root(x);	
	root_y = get_root(y);
	if(root_x != root_y){
		myunion(root_x, root_y);	
		flag = true;
	}
	return flag;
}
/* link_x and link_y will be constructed.*/
void make_network(int p){
	cout<<"input connectivity: ";
	cin>>p;
	for(int i=0;i<L; ++i){
		for(int j=0;j<L;++j){
		if(std_rand()<p){link_x[i][j]=1;}		
		if(std_rand()<p){link_y[i][j]=1;}		
	}}
}

void print_link(){
	cout<<"\nlink x direction."<<endl;
	for(int i=0;i<L; ++i){
	for(int j=0;j<L; ++j){
		if(link_x[i][j]==1){ cout<<"-"; }
		else if(link_x[i][j]==0){ cout<<" "; }
	}cout<<endl;
	}
	cout<<"\nlink y direction."<<endl;
	for(int i=0;i<L; ++i){
	for(int j=0;j<L; ++j){
		if(link_y[i][j]==1){ cout<<"|"; }
		else if(link_y[i][j]==0){ cout<<" "; }
	}cout<<endl;
	}
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

double E(int x, int y)
{
    // energy of bonds connected to site (x,y)
    // find the nearest neighbor sites with periodic boundary conditions
    int x_plus = x + 1, x_minus = x - 1, y_plus = y + 1, y_minus = y - 1;
    if (x_plus  ==  L)  x_plus  = 0;
    if (x_minus == -1)  x_minus = L - 1;
    if (y_plus  ==  L)  y_plus  = 0;
    if (y_minus == -1)  y_minus = L - 1;
    // contribution to energy from 4 bonds at site (x,y)
    double E = - S[x][y]*(
		   J_x[x_minus][y]*S[x_minus][y]+J_x[x][y]*S[x_plus][y] 
		 + J_y[x][y_minus]*S[x][y_minus]+J_y[x][y]*S[x][y_plus]   
		 );
    return E;
}

double E()
{
    double E_sum = 0;
    for (int x = 0; x < L; ++x)
        for (int y = 0; y < L; ++y)
            E_sum += E(x, y);
    return E_sum / 2;
}

double M()
{
    double M_sum = 0;
    for (int x = 0; x < L; ++x){
        for (int y = 0; y < L; ++y){
            M_sum += S[x][y];
	}}
    return M_sum ;
}

void print_S(){
	for(int i=0;i<L;++i){
		for(int j=0;j<L;++j){
			cout<<S[i][j]<<" ";
		}cout<<endl;
	}
}
//Total computational cost is 3*N = 3*L*L
void SW(double T){
/**** Link Cutting or Freezing Process. *******/	
	double p_f = 1.0 - exp(-2*J_0/T); // probability of frozen-bond. 
	int n_c=0; // number of label of claster.
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			set[i*L+j]=-1;
			link_x[i][j] = 0;	
			link_y[i][j] = 0;	
	}} // Initialize. This vector tells index of cluster.
	
	for(int x=0;x<L;x++){
		for(int y=0;y<L;y++){
			int x_plus=(x+1)%L; int y_plus=(y+1)%L;
			bool condtion_x = S[x][y]*S[x_plus][y] > 0;
			bool condtion_y = S[x][y]*S[x][y_plus] > 0;
			if( ( std_rand()<p_f ) && condtion_x ){ 
			link_x[x][y]=1; 
			myfind( x_plus *L+y, x*L+y);
			}
			if( ( std_rand()<p_f ) && condtion_y ){ 
			link_y[x][y]=1; 
			myfind( x*L+y, x*L+y_plus );
			}
		}}
/**** Update by using union-find algo. *******/	
	for(int k=0;k<N;++k){
		if(set[k]<0){ // is this order log N ??
			int state_of_cluster;
			double r=std_rand();
			if(r<=0.5){state_of_cluster=1;}
			if(r>0.5){state_of_cluster=-1;}
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
		}
		}//end of "set[k]<0" loop.
	}
}

int main(){
	int n_T = 100;	
	double T_max = 4.0,T_min = 1.0;
	double dT = (T_max - T_min)/ double(n_T);
	double e_step, e_av=0, e_sqd_av=0;      // average E per spin, squared average
	double m_step, m_av=0, m_sqd_av=0;      // average E per spin, squared average
	double c; double chi;
	int equilibration_steps = 3000;
	int production_steps = 20002;
	int np = 50;
	double accept_percent;
	double num_spins = double(L) * double(L);

    	cout << " Monte Carlo Simulation of the 2-dimensional Ising  Model\n"
		 << " ----------------------------------------------------\n"
		 << " " << L*L << " spins on " << L << "x" << L << " lattice," << endl;

	cout << "T\t"<< "<e>\t"<< "<e^2>\t"<< "<m>\t"<< "<m^2>"<< "c\t"<< "chi"  << endl;
	 for(int k = 0; k<n_T;k++){
	    	T = T_min + k*dT;	
		initialize();
		init_S_neel();	// set spins, which are all up.
		init_J_pure_ferro();	// set coupling constants for all 1.0
		
		for(int step = 0; step < equilibration_steps; ++step) {
			SW(T);		
		}
		e_av=0, e_sqd_av=0;
		m_av=0, m_sqd_av=0;
		for(int step=0; step < production_steps; ++step){
			if(step%np == 0){
				SW(T);		
				m_step = M() / num_spins;
				e_step = E() / num_spins;
				e_av += e_step;
				e_sqd_av += e_step * e_step;
				m_av += m_step;
				m_sqd_av += m_step * m_step;
				}
		}
		double per_step= double(np) / double(production_steps);
		e_av *= per_step;
		e_sqd_av *= per_step;
		m_av *= per_step;
		m_sqd_av *= per_step;
		c = (e_sqd_av - e_av * e_av) / (T * T);
		chi = (m_sqd_av - m_av * m_av) / T;
		cout << T <<'\t'<< e_av << '\t'<< e_sqd_av <<'\t'
			<< m_av << '\t'<< m_sqd_av << '\t'<< c <<'\t' << chi << endl;
	}
	return 0;
}
